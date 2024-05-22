import numpy as np
from networkx import DiGraph
import networkx as nx
from itertools import combinations

from subgraphs.sub_graphs_utils import get_adjacency_matrix
from utils.types import LargeSubGraphSearchResult


class SingleInputModule:
    """
    SIM - Single Input Module. my variation of induced SIM detector.
    paper:    Shai S. Shen-Orr1, Ron Milo2, Shmoolik Mangan1 & Uri Alon1
              "Network motifs in the transcriptional regulation  network of Escherichia coli"
    """

    def __init__(self, network: DiGraph):
        self.network = network

        self.fsl = {}  # frequent sub graph list - value is the frequency
        self.fsl_fully_mapped = {}  # same fsl, the value is the list of sub graphs
        self.adj_mats = {}

    def search_sub_graphs(self, min_control_size: int) -> LargeSubGraphSearchResult:

        _, max_out_degree = max(self.network.out_degree, key=lambda x: x[1])
        # we don't want the keys to mix with regular motif ids
        self.fsl = {f'SIM_{i}': 0 for i in range(min_control_size, max_out_degree + 1)}
        self.fsl_fully_mapped = {f'SIM_{i}': [] for i in range(min_control_size, max_out_degree + 1)}
        self.adj_mats = {f'SIM_{i}': np.array([]) for i in range(min_control_size, max_out_degree + 1)}

        for control_size in range(min_control_size, max_out_degree + 1):
            sim_key = f'SIM_{control_size}'
            for input_node in list(self.network.nodes):
                neighbors = list(self.network.neighbors(input_node))
                for controlled in list(combinations(neighbors, control_size)):
                    sub_graph = list(controlled) + [input_node]
                    induced_sim = nx.induced_subgraph(self.network, sub_graph)
                    if len(induced_sim.edges) == control_size:
                        self.fsl[sim_key] += 1
                        self.fsl_fully_mapped[sim_key].append(tuple(list(induced_sim.edges)))

            if self.fsl[sim_key]:
                sub_graph = list(self.fsl_fully_mapped[sim_key][0])
                adj_mat = get_adjacency_matrix(nx.DiGraph(sub_graph))
                self.adj_mats[sim_key] = adj_mat

        return LargeSubGraphSearchResult(fsl=self.fsl,
                                         fsl_fully_mapped=self.fsl_fully_mapped,
                                         adj_mat=self.adj_mats)
