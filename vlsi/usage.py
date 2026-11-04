from microfluidic_pr_bridge import MicrofluidicPRBridge
from MMFT import simulateMMFTWithNodeData

# 1. Initialize the Place and Route bridge with your droplet circuit JSON specification
pr_bridge = MicrofluidicPRBridge(json_spec_path="exemplo_3.json")

# 2. Run Simulated Annealing for optimal cell placement
print("Running Simulated Annealing Placement...")
optimal_positions = pr_bridge.run_placement(temperatura_inicial=400.0, iteracoes=100)

print(f"Placed {len(optimal_positions)} microfluidic nodes.")

# 3. Run Ant Colony Optimization (ACO) for global channel routing
print("Running ACO Global Router...")
router = pr_bridge.run_routing(router_type="aco", n_iteracoes=30)

metrics = router.metricas()
print(f"Routing Complete -> Wirelength: {metrics['wirelength']}, Overflow: {metrics['overflow']}")

# 4. Export the P&R optimized physical layout directly into MMFT simulation format
mmft_input_data = pr_bridge.export_to_mmft_simulation_format()

# 5. Simulate using the Munich Microfluidic Toolkit simulator
print("Executing MMFT Simulator with optimized P&R layout...")
simulation_results = simulateMMFTWithNodeData(mmft_input_data)

print("Simulation finished successfully!")



from cts_dme_placement import CTSDMEMergeClient, CTSDMEPlacer
from microfluidic_pr_bridge import MicrofluidicPRBridge

# 1. Load microfluidic project specification
spec_path = "exemplo_3.json"

# 2. Run CTS + DME to determine optimal synchronized pump positions on the chip perimeter
placer = CTSDMEPlacer(json_spec_path=spec_path if 'json_spec_path' in globals() else None, chip_width=50.0, chip_height=50.0)
# Alternatively, pass loaded dict data directly:
import json
with open("exemplo_3.json", "r") as f:
    data = json.load(f)

cts_placer = CTSDMEPlacer(data, chip_width=40.0, chip_height=40.0)
pump_coordinates = cts_placer.run_cts_dme()

print("Synchronized Pump & Junction Coordinates from CTS-DME:")
for node_id, coords in pump_coordinates.items():
    print(f"  Node/Pump {node_id}: {coords}")

# 3. Feed coordinates and optimized layout into the P&R Bridge
pr_bridge = MicrofluidicPRBridge(data_dict=data)
pr_bridge.run_placement() # Refines overall placement
router = pr_bridge.run_routing(router_type="aco") # Routes channels avoiding congestion

# 4. Export to MMFT simulator
mmft_input = pr_bridge.export_to_mmft_simulation_format()