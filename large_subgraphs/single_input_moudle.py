import numpy as np
from networkx import DiGraph
import networkx as nx
from itertools import combinations
from utils.types import LargeSubGraphSearchResult


class SingleInputModule:
    """
    SIM - Single Input Module. my variation of induced SIM detector.
    paper:    Shai S. Shen-Orr1, Ron Milo2, Shmoolik Mangan1 & Uri Alon1
              "Network motifs in the transcriptional regulation  network of Escherichia coli"
    """

    def __init__(self, network: DiGraph):
        self.network = network
        s, t = list(network.edges)[0]
        self.use_polarity = 'polarity' in network[s][t] and network[s][t]['polarity'] is not None

        self.fsl = {}  # frequent sub graph list - value is the frequency
        self.fsl_fully_mapped = {}  # same fsl, the value is the list of sub graphs
        self.adj_mats = {}

    @staticmethod
    def __set_adj_mat(control_size: int) -> np.ndarray:
        adj_mat = np.zeros((control_size + 1, control_size + 1))
        adj_mat[0][1:control_size + 1] = 1
        return adj_mat

    def search_sub_graphs(self, min_control_size: int, max_control_size: int) -> LargeSubGraphSearchResult:
        _, max_out_degree = max(self.network.out_degree, key=lambda x: x[1])

        if max_control_size is None:
            max_control_size = min(max_out_degree, max_out_degree)

        # we don't want the keys to mix with regular motif ids
        self.fsl = {f'SIM_{i}': 0 for i in range(min_control_size, max_control_size + 1)}
        self.fsl_fully_mapped = {f'SIM_{i}': [] for i in range(min_control_size, max_control_size + 1)}
        self.adj_mats = {f'SIM_{i}': self.__set_adj_mat(i) for i in range(min_control_size, max_control_size + 1)}

        for control_size in range(min_control_size, max_control_size + 1):
            sim_key = f'SIM_{control_size}'
            for input_node in list(self.network.nodes):
                neighbors = list(self.network.neighbors(input_node))
                if input_node in neighbors:
                    neighbors.remove(input_node)
                for controlled in list(combinations(neighbors, control_size)):
                    sub_graph = list(controlled) + [input_node]
                    induced_sim = nx.induced_subgraph(self.network, sub_graph)
                    if len(induced_sim.edges) != control_size:
                        continue

                    self.fsl[sim_key] += 1

                    if not self.use_polarity:
                        self.fsl_fully_mapped[sim_key].append(tuple(list(induced_sim.edges)))
                    else:
                        polarities = nx.get_edge_attributes(induced_sim, 'polarity')
                        pol_edges = tuple([(*e, {'polarity': polarities[e]}) for e in list(induced_sim.edges)])
                        self.fsl_fully_mapped[sim_key].append(pol_edges)

        return LargeSubGraphSearchResult(fsl=self.fsl,
                                         fsl_fully_mapped=self.fsl_fully_mapped,
                                         adj_mat=self.adj_mats)
