# -*- coding: utf-8 -*-
"""
cts_dme_placement.py
====================
Clock Tree Synthesis (CTS) & Deferred Merge Embedding (DME) for 
Microfluidic Pump Placement and Droplet Synchronization.
"""

from __future__ import annotations
import math
import numpy as np
from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Optional


@dataclass
class Segment2D:
    """Represents a DME locus (a line segment in 2D space) for a merge node."""
    p1: np.ndarray  # [x, y]
    p2: np.ndarray  # [x, y]

    def midpoint(self) -> np.ndarray:
        return (self.p1 + self.p2) / 2.0


@dataclass
class TreeNode:
    id: str
    is_leaf: bool = False
    left: Optional[TreeNode] = None
    right: Optional[TreeNode] = None
    parent: Optional[TreeNode] = None
    target_time: float = 0.0
    position: Optional[np.ndarray] = None
    locus: Optional[Segment2D] = None
    radius: float = 0.0  # Allowed skew/delay buffer radius


class CTSDMEPlacer:
    """
    Implements Deferred Merge Embedding (DME) combined with Clock Tree Synthesis 
    to position microfluidic merge junctions and anchor input pumps.
    """

    def __init__(
        self,
        mix_data: dict,
        chip_width: float = 100.0,
        chip_height: float = 100.0,
        velocity: float = 1.0  # droplet propagation speed in units/sec
    ):
        self.mix_data = mix_data
        self.chip_width = chip_width
        self.chip_height = chip_height
        self.velocity = velocity
        self.nodes: Dict[str, TreeNode] = {}
        self.root: Optional[TreeNode] = None

    def build_mixing_tree(self):
        """Constructs the binary mixing tree from MMFT/CRN specification."""
        initials = self.mix_data.get("goticulasIniciais", [])
        mixes = self.mix_data.get("goticulasMisturadas", [])

        # Create leaf nodes for initial droplets (Pumps)
        for d in initials:
            did = str(d.get("id"))
            self.nodes[did] = TreeNode(id=did, is_leaf=True, target_time=0.0)

        # Create internal nodes for mixed droplets
        for m in mixes:
            mid = str(m.get("id"))
            parents = [str(p) for p in m.get("goticulasMisturadas", [])]
            mix_times = m.get("tempoMistura", [0.0] * len(parents))
            cumulative_time = sum(mix_times)

            node = TreeNode(id=mid, is_leaf=False, target_time=cumulative_time)
            self.nodes[mid] = node

            # Link parents to the binary merge node
            if len(parents) >= 2:
                left_id, right_id = parents[0], parents[1]
                if left_id in self.nodes and right_id in self.nodes:
                    node.left = self.nodes[left_id]
                    node.right = self.nodes[right_id]
                    self.nodes[left_id].parent = node
                    self.nodes[right_id].parent = node

            self.root = node  # The last created mix is the root/sink

    def _compute_locus_bottom_up(self, node: Optional[TreeNode]):
        """Phase 1: Bottom-Up Locus Generation."""
        if node is None or node.is_leaf:
            if node and node.is_leaf:
                # Leaf locus is a point on the chip boundary (assigned later or initialized near center)
                pos = np.array([np.random.uniform(10, self.chip_width - 10), 0.0])
                node.position = pos
                node.locus = Segment2D(pos, pos)
            return

        self._compute_locus_bottom_up(node.left)
        self._compute_locus_bottom_up(node.right)

        if node.left and node.right:
            l1 = node.left.locus.midpoint()
            l2 = node.right.locus.midpoint()
            
            # The merge locus is a segment connecting child loci centers adjusted for arrival time delta
            mid = (l1 + l2) / 2.0
            direction = l2 - l1
            length = np.linalg.norm(direction)
            
            if length > 1e-5:
                direction = direction / length
                perp = np.array([-direction[1], direction[0]])
            else:
                perp = np.array([1.0, 0.0])

            # Allow a small segment length based on timing skew flexibility
            time_delta = abs(node.left.target_time - node.right.target_time)
            segment_radius = max(5.0, time_delta * self.velocity * 0.5)
            
            p1 = mid - perp * segment_radius
            p2 = mid + perp * segment_radius
            node.locus = Segment2D(p1, p2)
            node.position = mid

    def _embed_top_down(self, node: Optional[TreeNode], inherited_pos: Optional[np.ndarray] = None):
        """Phase 2: Top-Down Coordinate Embedding and Pump Anchoring."""
        if node is None:
            return

        if inherited_pos is not None:
            # Project inherited parent position onto current node's locus segment
            p1, p2 = node.locus.p1, node.locus.p2
            vec = p2 - p1
            len_sq = np.dot(vec, vec)
            if len_sq < 1e-5:
                node.position = p1
            else:
                t = np.dot(inherited_pos - p1, vec) / len_sq
                t = np.clip(t, 0.0, 1.0)
                node.position = p1 + t * vec
        elif node.locus:
            node.position = node.locus.midpoint()

        # If it's a leaf node, snap/anchor it to the nearest chip boundary edge (Pump position)
        if node.is_leaf and node.position is not None:
            x, y = node.position
            # Find closest border (Left, Right, Top, Bottom)
            distances = [x, self.chip_width - x, y, self.chip_height - y]
            min_idx = np.argmin(distances)
            
            if min_idx == 0:   node.position = np.array([0.0, np.clip(y, 0, self.chip_height)])
            elif min_idx == 1: node.position = np.array([self.chip_width, np.clip(y, 0, self.chip_height)])
            elif min_idx == 2: node.position = np.array([np.clip(x, 0, self.chip_width), 0.0])
            elif min_idx == 3: node.position = np.array([np.clip(x, 0, self.chip_width), self.chip_height])

        if node.left:
            self._embed_top_down(node.left, node.position)
        if node.right:
            self._embed_top_down(node.right, node.position)

    def run_cts_dme(self) -> Dict[str, Tuple[float, float]]:
        """
        Executes the CTS + DME algorithm and returns precise pump and junction coordinates.
        """
        self.build_mixing_tree()
        if not self.root:
            raise ValueError("Invalid mixing tree: No root node found.")
        
        self._compute_locus_bottom_up(self.root)
        self._embed_top_down(self.root)

        # Collect final coordinates for all nodes (pumps & junctions)
        coordinates = {}
        for nid, node in self.nodes.items():
            if node.position is not None:
                coordinates[nid] = (float(node.position[0]), float(node.position[1]))
        return coordinates