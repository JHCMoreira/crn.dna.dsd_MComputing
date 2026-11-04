"""Netlist-to-chip CAD compiler.

The compiler intentionally has no dependency on the MMFT engine.  It produces
the v2 project format consumed by the existing simulator, which makes the
placement and routing result inspectable before a simulation is submitted.
"""

from __future__ import annotations

from dataclasses import dataclass
import heapq
import math
from typing import Any


@dataclass(frozen=True)
class Pin:
    component: str
    x: float | None = None
    y: float | None = None
    role: str = "junction"


@dataclass(frozen=True)
class Net:
    id: str
    pins: tuple[Pin, ...]
    target_time: float = 0.0


def _id(value: Any, prefix: str, index: int) -> str:
    text = str(value if value is not None else f"{prefix}{index}").strip()
    if not text:
        raise ValueError(f"{prefix} identifiers cannot be empty")
    return text


def _normalise(data: dict[str, Any]) -> tuple[list[Pin], list[Net], dict[str, Any]]:
    raw_components = data.get("junctions", data.get("components", data.get("nodes", [])))
    if not isinstance(raw_components, list):
        raise ValueError("junctions/components/nodes must be a list")
    components: dict[str, Pin] = {}
    for i, item in enumerate(raw_components):
        if isinstance(item, str):
            item = {"id": item}
        if not isinstance(item, dict):
            raise ValueError("each junction must be an object or identifier")
        cid = _id(item.get("id", item.get("name")), "j", i)
        if cid in components:
            raise ValueError(f"duplicate junction: {cid}")
        role = str(item.get("role", item.get("type", "junction"))).lower()
        components[cid] = Pin(cid, item.get("x"), item.get("y"), role)

    raw_nets = data.get("nets", data.get("netlist", []))
    if not isinstance(raw_nets, list):
        raise ValueError("nets/netlist must be a list")
    nets: list[Net] = []
    for i, item in enumerate(raw_nets):
        if not isinstance(item, dict):
            raise ValueError("each net must be an object")
        nid = _id(item.get("id", item.get("name")), "net", i)
        raw_pins = item.get("pins", item.get("junctions", item.get("components", [])))
        if not isinstance(raw_pins, list) or len(raw_pins) < 2:
            raise ValueError(f"{nid}: a net requires at least two pins")
        pins: list[Pin] = []
        for pin in raw_pins:
            cid = str(pin.get("id", pin.get("name"))) if isinstance(pin, dict) else str(pin)
            if cid not in components:
                components[cid] = Pin(cid)
            pins.append(components[cid])
        nets.append(Net(nid, tuple(pins), float(item.get("target_time", item.get("time", 0.0)))))

    # A channel list is a convenient netlist shorthand.  Each channel becomes
    # a two-pin net unless an equivalent net was supplied explicitly.
    for i, channel in enumerate(data.get("channels", [])):
        if not isinstance(channel, dict):
            raise ValueError("each channel must be an object")
        source = str(channel.get("source", channel.get("from")))
        target = str(channel.get("target", channel.get("to")))
        if source not in components:
            components[source] = Pin(source)
        if target not in components:
            components[target] = Pin(target)
        nets.append(Net(_id(channel.get("id"), "channel", i),
                        (components[source], components[target]),
                        float(channel.get("target_time", channel.get("time", 0.0)))))
    if not nets:
        raise ValueError("netlist requires nets or channels")
    return list(components.values()), nets, data


def _placement(components: list[Pin], nets: list[Net], grid_slack: float) -> dict[str, tuple[int, int]]:
    count = len(components)
    columns = max(2, math.ceil(math.sqrt(count * max(1.0, grid_slack))))
    positions: dict[str, tuple[int, int]] = {}
    anchored = [c for c in components if c.x is not None and c.y is not None]
    for component in anchored:
        positions[component.component] = (round(float(component.x)), round(float(component.y)))
    unplaced = [c for c in components if c.component not in positions]
    # Place high-connectivity junctions first near the centroid of their pins.
    degree = {c.component: sum(c.component in (p.component for p in n.pins) for n in nets) for c in components}
    unplaced.sort(key=lambda c: (-degree[c.component], c.component))
    occupied = set(positions.values())
    for component in unplaced:
        neighbours = [positions[p.component] for n in nets if component.component in (p.component for p in n.pins)
                      for p in n.pins if p.component != component.component and p.component in positions]
        if neighbours:
            center = (round(sum(p[0] for p in neighbours) / len(neighbours)),
                      round(sum(p[1] for p in neighbours) / len(neighbours)))
            candidates = sorted(((abs(x - center[0]) + abs(y - center[1]), x, y)
                                 for y in range(columns + 1) for x in range(columns + 1)
                                 if (x, y) not in occupied))
            _, x, y = candidates[0]
        else:
            slot = len(occupied)
            x, y = slot % columns, slot // columns
            while (x, y) in occupied:
                slot += 1
                x, y = slot % columns, slot // columns
        positions[component.component] = (x, y)
        occupied.add((x, y))
    return positions


def _route(start: tuple[int, int], goal: tuple[int, int], blocked: set[tuple[int, int]],
           capacity: dict[tuple[tuple[int, int], tuple[int, int]], int],
           width: int, height: int) -> list[tuple[int, int]]:
    if start == goal:
        return [start]
    queue = [(0, 0, start)]
    previous: dict[tuple[int, int], tuple[int, int] | None] = {start: None}
    distance = {start: 0}
    while queue:
        _, _, current = heapq.heappop(queue)
        if current == goal:
            break
        for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            nxt = (current[0] + dx, current[1] + dy)
            if not (0 <= nxt[0] <= width and 0 <= nxt[1] <= height) or nxt in blocked and nxt != goal:
                continue
            edge = tuple(sorted((current, nxt)))
            cost = 1 + capacity.get(edge, 0) * 4
            candidate = distance[current] + cost
            if candidate < distance.get(nxt, math.inf):
                distance[nxt] = candidate
                previous[nxt] = current
                heapq.heappush(queue, (candidate + abs(nxt[0] - goal[0]) + abs(nxt[1] - goal[1]),
                                       len(previous), nxt))
    if goal not in previous:
        raise ValueError(f"unable to route {start} to {goal}")
    path = []
    cursor: tuple[int, int] | None = goal
    while cursor is not None:
        path.append(cursor)
        cursor = previous[cursor]
    return list(reversed(path))


def compile_netlist(data: dict[str, Any]) -> dict[str, Any]:
    """Compile a netlist into a validated v2 project plus CAD diagnostics."""
    if not isinstance(data, dict):
        raise ValueError("CAD input must be an object")
    components, nets, source = _normalise(data)
    grid_slack = float(source.get("grid_slack", 1.3))
    grid_pitch = float(source.get("grid_pitch", 1e-3))
    if not math.isfinite(grid_pitch) or grid_pitch <= 0:
        raise ValueError("grid_pitch must be positive")
    capacity_limit = max(1, int(source.get("routing_capacity", 3)))
    positions = _placement(components, nets, grid_slack)
    max_x = max(x for x, _ in positions.values())
    max_y = max(y for _, y in positions.values())
    blocked = set(positions.values())
    usage: dict[tuple[tuple[int, int], tuple[int, int]], int] = {}
    routes: list[dict[str, Any]] = []
    for net in nets:
        # CTS/DME equivalent: grow a shortest-path tree from the median pin,
        # attaching remaining pins in nearest-first order.
        terminals = [positions[p.component] for p in net.pins]
        tree = [terminals[0]]
        remaining = terminals[1:]
        segments: list[list[tuple[int, int]]] = []
        while remaining:
            goal = min(remaining, key=lambda point: min(abs(point[0] - t[0]) + abs(point[1] - t[1]) for t in tree))
            start = min(tree, key=lambda point: abs(point[0] - goal[0]) + abs(point[1] - goal[1]))
            path = _route(start, goal, blocked - {start, goal}, usage, max_x, max_y)
            for a, b in zip(path, path[1:]):
                edge = tuple(sorted((a, b)))
                usage[edge] = usage.get(edge, 0) + 1
                if usage[edge] > capacity_limit:
                    raise ValueError(f"routing capacity exceeded on net {net.id}")
            tree.extend(path)
            remaining.remove(goal)
            segments.append(path)
        routes.append({"id": net.id, "pins": [p.component for p in net.pins], "segments": segments,
                       "target_time": net.target_time})

    node_ids = {c.component: f"j{i}" for i, c in enumerate(components)}
    coordinate_ids = {positions[c.component]: node_ids[c.component] for c in components}
    route_coordinates = {point for route in routes for segment in route["segments"] for point in segment}
    for point in sorted(route_coordinates):
        coordinate_ids.setdefault(point, f"r{len(coordinate_ids)}")
    component_by_id = {node_ids[c.component]: c for c in components}
    nodes = []
    for point, node_id in sorted(coordinate_ids.items(), key=lambda item: item[1]):
        component = component_by_id.get(node_id)
        nodes.append({"id": node_id, "name": component.component if component else node_id,
                      "x": point[0] * grid_pitch, "y": point[1] * grid_pitch,
                      "sink": bool(component and component.role in {"sink", "output"}),
                      "ground": bool(component and component.role in {"sink", "output"})})
    channels: list[dict[str, Any]] = []
    route_channel_ids: dict[str, list[str]] = {}
    for route in routes:
        route_channel_ids[route["id"]] = []
        for segment_index, segment in enumerate(route["segments"]):
            for part, (a, b) in enumerate(zip(segment, segment[1:])):
                a_id = coordinate_ids[a]
                b_id = coordinate_ids[b]
                channel_id = f"{route['id']}:c{segment_index}:{part}"
                channels.append({"id": channel_id,
                                 "source": a_id, "target": b_id,
                                 "width": float(source.get("channel_width", 50e-6)),
                                 "height": float(source.get("channel_height", 30e-6))})
                route_channel_ids[route["id"]].append(channel_id)
    sinks = [c for c in components if c.role in {"sink", "output"}]
    if not sinks:
        sinks = [max(components, key=lambda c: positions[c.component][1])]
    sink_id = node_ids[sinks[0].component]
    pumps = [{"id": f"pump:{p.component}", "source": node_ids[p.component], "target": sink_id,
              "kind": "flow", "value": float(source.get("pump_flow", 1.5e-16))}
             for p in components if p.role in {"source", "input"} and node_ids[p.component] != sink_id]
    if not pumps:
        inlet = next((c for c in components if node_ids[c.component] != sink_id), None)
        if inlet is not None:
            pumps = [{"id": "pump:inlet", "source": node_ids[inlet.component], "target": sink_id,
                      "kind": "flow", "value": float(source.get("pump_flow", 1.5e-16))}]
    raw_droplets = source.get("droplets", source.get("injections", []))
    if not isinstance(raw_droplets, list):
        raise ValueError("droplets/injections must be a list")
    species_keys = sorted({str(key) for droplet in raw_droplets if isinstance(droplet, dict)
                           for key in droplet.get("concentrations", {})})
    species_ids = {key: f"s{i}" for i, key in enumerate(species_keys)}
    species = [{"id": species_ids[key], "name": f"species_{i}"} for i, key in enumerate(species_keys)]
    injections = []
    for i, droplet in enumerate(raw_droplets):
        if not isinstance(droplet, dict):
            raise ValueError("each droplet must be an object")
        route_id = str(droplet.get("net", droplet.get("channel", "")))
        channel_id = droplet.get("channel")
        if channel_id not in {channel["id"] for channel in channels}:
            channel_id = route_channel_ids.get(route_id, [None])[0]
        if channel_id is None:
            raise ValueError(f"droplet {i}: channel or net does not identify a routed channel")
        concentrations = {species_ids[str(key)]: float(value)
                          for key, value in droplet.get("concentrations", {}).items()}
        injections.append({"id": _id(droplet.get("id"), "droplet", i),
                            "name": str(droplet.get("name", droplet.get("id", f"droplet-{i}"))),
                            "channel": channel_id, "time": float(droplet.get("time", 0.1)),
                            "position": float(droplet.get("position", 0.25)),
                            "volume": float(droplet.get("volume", 2.25e-13)),
                            "concentrations": concentrations})
    project = {"schema_version": 2, "name": source.get("name", "CAD circuit"), "nodes": nodes,
               "channels": channels, "pumps": pumps, "species": species, "injections": injections,
               "provenance": {"profile": "place-route-cts-dme",
               "netlist": source}}
    return {"project": project, "placement": {c.component: {"x": positions[c.component][0],
            "y": positions[c.component][1]} for c in components}, "routes": routes,
            "cts": {"algorithm": "CTS/DME", "sink": sinks[0].component, "capacity": capacity_limit}}
