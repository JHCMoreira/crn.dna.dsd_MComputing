import pytest

from backend.cad import compile_netlist
from backend.models import Project


def test_compile_netlist_places_routes_and_builds_project():
    result = compile_netlist({
        "name": "mix",
        "junctions": [
            {"id": "a", "role": "source"},
            {"id": "b", "role": "source"},
            {"id": "mix", "role": "sink"},
        ],
        "nets": [{"id": "mix-net", "pins": ["a", "b", "mix"], "target_time": 2.5}],
    })

    project = Project.model_validate(result["project"])
    assert set(result["placement"]) == {"a", "b", "mix"}
    assert result["routes"][0]["target_time"] == 2.5
    assert project.channels
    assert len(project.pumps) == 2
    assert result["cts"]["algorithm"] == "CTS/DME"


def test_channels_are_accepted_as_netlist_shorthand():
    result = compile_netlist({
        "junctions": [{"id": "in", "role": "source"}, {"id": "out", "role": "sink"}],
        "channels": [{"id": "feed", "source": "in", "target": "out"}],
    })
    assert result["routes"][0]["id"] == "feed"
    Project.model_validate(result["project"])


def test_droplets_are_mapped_to_routed_channels():
    result = compile_netlist({
        "junctions": [{"id": "in", "role": "source"}, {"id": "out", "role": "sink"}],
        "nets": [{"id": "feed", "pins": ["in", "out"]}],
        "droplets": [{"id": "sample", "net": "feed", "concentrations": {"DNA": 2.0}}],
    })
    injection = result["project"]["injections"][0]
    assert injection["channel"].startswith("feed:")
    assert injection["concentrations"] == {"s0": 2.0}


def test_invalid_netlist_is_rejected():
    with pytest.raises(ValueError, match="at least two pins"):
        compile_netlist({"junctions": [{"id": "a"}], "nets": [{"id": "bad", "pins": ["a"]}]})
