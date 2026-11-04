# -*- coding: utf-8 -*-
"""
microfluidic_pr_bridge.py
=========================
Library that bridges VLSI Place-and-Route algorithms (Simulated Annealing 
Placement + ACO Global Routing) with the Munich Microfluidic Toolkit (MMFT) 
droplet simulator.
"""

from __future__ import annotations
import json
import math
from typing import Dict, List, Tuple, Optional
import numpy as np
import pandas as pd

# Import algorithms from annexed modules (assumed available in PYTHONPATH)
from sa_vlsi_placement import (
    SAPlacement,
    carregar_componentes,
    carregar_netlist,
    construir_redes,
    dimensoes_grade_sugerida
)
from aco_vlsi_routing import (
    GradeRoteamento,
    construir_pinos_das_redes,
    RoteadorACO,
    RoteadorMaze,
    RoteadorPadrao
)


class MicrofluidicPRBridge:
    """
    Orchestrator for Place and Route on Microfluidic Droplet Biochips,
    translating abstract logical CRNs/mixings into physical MMFT simulation inputs.
    """

    def __init__(
        self,
        json_spec_path: Optional[str] = None,
        data_dict: Optional[dict] = None,
        grid_slack: float = 1.3,
        routing_capacity: float = 3.0
    ):
        if data_dict is not None:
            self.data = data_dict
        elif json_spec_path is not None:
            with open(json_spec_path, 'r', encoding='utf-8') as f:
                self.data = json.load(f)
        else:
            raise ValueError("Provide either json_spec_path or data_dict.")

        self.grid_slack = grid_slack
        self.routing_capacity = routing_capacity
        
        self.df_componentes: Optional[pd.DataFrame] = None
        self.df_netlist: Optional[pd.DataFrame] = None
        self.redes: Dict[str, List[int]] = {}
        self.n_linhas = 0
        self.n_colunas = 0
        self.posicoes: Dict[int, Tuple[int, int]] = {}
        self.grade: Optional[GradeRoteamento] = None
        self.router_result = None

    def parse_microfluidic_to_netlist(self) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Converts MMFT-style droplet specification ('goticulasIniciais' and 
        'goticulasMisturadas') into standard component and netlist DataFrames
        compatible with the VLSI P&R engines.
        """
        component_names = set()
        initial_droplets = self.data.get("goticulasIniciais", [])
        mixed_droplets = self.data.get("goticulasMisturadas", [])

        # Collect all component/droplet identifiers
        for d in initial_droplets:
            component_names.add(str(d.get("id")))
        for m in mixed_droplets:
            component_names.add(str(m.get("id")))
            for comp_id in m.get("goticulasMisturadas", []):
                component_names.add(str(comp_id))

        comp_list = sorted(list(component_names))
        df_comp = pd.DataFrame({
            "Componente": comp_list,
            "Largura": [1.0] * len(comp_list),
            "Altura": [1.0] * len(comp_list)
        })

        # Build netlist from mixing dependencies
        netlist_rows = []
        for m in mixed_droplets:
            net_id = f"Net_{m.get('id')}"
            parent_droplets = m.get("goticulasMisturadas", [])
            # Include the mixed product and its parents in the net
            net_members = set([str(m.get("id"))] + [str(p) for p in parent_droplets])
            for member in net_members:
                netlist_rows.append({"Net": net_id, "Componente": member})

        df_net = pd.DataFrame(netlist_rows)
        
        self.df_componentes = df_comp
        self.df_netlist = df_net
        return df_comp, df_net

    def run_placement(self, temperatura_inicial: float = 500.0, iteracoes: int = 150) -> Dict[int, Tuple[int, int]]:
        """
        Executes Simulated Annealing (SA) to position microfluidic nodes on a grid.
        """
        if self.df_componentes is None or self.df_netlist is None:
            self.parse_microfluidic_to_netlist()

        n_comps = len(self.df_componentes)
        self.n_linhas, self.n_colunas = dimensoes_grade_sugerida(n_comps, folga=self.grid_slack)
        self.redes = construir_redes(self.df_componentes, self.df_netlist)

        sa = SAPlacement(
            n_componentes=n_comps,
            redes=self.redes,
            n_linhas=self.n_linhas,
            n_colunas=self.n_colunas,
            temperatura_inicial=temperatura_inicial,
            iteracoes_por_temperatura=iteracoes,
            semente=42
        )
        historico = sa.executar()
        self.posicoes = historico.melhor_posicoes_global
        return self.posicoes

    def run_routing(self, router_type: str = "aco", n_iteracoes: int = 50) -> RoteadorGlobal:
        """
        Executes global routing (ACO, Maze, or L-Pattern) over the placed microfluidic grid.
        """
        if not self.posicoes:
            self.run_placement()

        self.grade = GradeRoteamento(self.n_linhas, self.n_colunas, capacidade=self.routing_capacity)
        pinos_redes = construir_pinos_das_redes(self.redes, self.posicoes, self.grade)

        if router_type.lower() == "aco":
            router = RoteadorACO(self.grade, pinos_redes, n_formigas=6, semete=42)
        elif router_type.lower() == "maze":
            router = RoteadorMaze(self.grade, pinos_redes, semente=42)
        else:
            router = RoteadorPadrao(self.grade, pinos_redes, semente=42)

        self.router_result = router
        router.executar(n_iteracoes=n_iteracoes)
        return router

    def export_to_mmft_simulation_format(self, base_scale_m: float = 1e-5) -> dict:
        """
        Translates the P&R grid placement and routing results into the exact JSON
        schema expected by MMFT's `simulateMMFTWithNodeData`.
        """
        if not self.posicoes or self.grade is None:
            raise RuntimeError("Run placement and routing before exporting to MMFT format.")

        nodes_out = []
        # Add a designated sink node
        sink_node_id = "Sink"
        sink_grid_pos = (self.n_linhas - 1, self.n_colunas - 1)

        # Build nodes from placement positions
        for idx, (linha, coluna) in self.posicoes.items():
            comp_name = str(self.df_componentes.loc[idx, "Componente"])
            is_sink_node = (linha == sink_grid_pos[0] and coluna == sink_grid_pos[1])
            
            nodes_out.append({
                "name": comp_name,
                "x": float(coluna),
                "y": float(linha),
                "isSink": is_sink_node,
                "ground": is_sink_node
            })

        # Ensure Sink node exists in nodes list
        if not any(n["name"] == sink_node_id for n in nodes_out):
            nodes_out.append({
                "name": sink_node_id,
                "x": float(sink_grid_pos[1]),
                "y": float(sink_grid_pos[0]),
                "isSink": True,
                "ground": True
            })

        # Build edges (channels) from routing results
        edges_out = []
        name_to_node = {n["name"]: n for n in nodes_out}

        # Map grid node index to name
        node_idx_to_name = {}
        for idx, (linha, coluna) in self.posicoes.items():
            node_idx_to_name[self.grade.no(linha, coluna)] = str(self.df_componentes.loc[idx, "Componente"])

        # Add routing paths as channels
        edge_counter = 0
        for i, net_name in enumerate(self.router_result.nomes):
            route = self.router_result.rotas[i]
            if route is None:
                continue
            for edge_id in route:
                u, v = self.grade.arestas[edge_id]
                u_name = node_idx_to_name.get(u, f"Node_{u}")
                v_name = node_idx_to_name.get(v, f"Node_{v}")

                if u_name not in name_to_node or v_name not in name_to_node:
                    continue

                edges_out.append({
                    "from": {"name": u_name},
                    "to": {"name": v_name},
                    "height": 3.0  # Normalized channel height parameter
                })
                edge_counter += 1

        # Fallback: if router didn't populate edges, connect components to sink
        if not edges_out:
            for node in nodes_out:
                if not node["isSink"]:
                    edges_out.append({
                        "from": {"name": node["name"]},
                        "to": {"name": sink_node_id},
                        "height": 3.0
                    })

        # Format initial droplets / pumps from input data
        initial_droplets_out = []
        for goticula in self.data.get("goticulasIniciais", []):
            initial_droplets_out.append({
                "id": goticula.get("id"),
                "nome": goticula.get("nome", goticula.get("id")),
                "pump": str(goticula.get("pump", goticula.get("id")))
            })

        mmft_payload = {
            "nodes": nodes_out,
            "edges": edges_out,
            "goticulasIniciais": initial_droplets_out,
            "goticulasMisturadas": self.data.get("goticulasMisturadas", [])
        }

        return mmft_payload