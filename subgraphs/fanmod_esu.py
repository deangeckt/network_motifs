import random

from networkx import DiGraph

from subgraphs.sub_graphs_abc import SubGraphsABC
import networkx as nx

from utils.sub_graphs import graph_to_hashed_graph
from collections import defaultdict

from utils.types import SubGraphSearchResult


class FanmodESU(SubGraphsABC):
    """
    FanMOD / ESU
    S. Wernicke, “Efficient detection of network motifs,” 2006
    pseudocode: Ribeiro, Pedro and Silva, Fernando and Kaiser, Marcus: "Strategies for Network Motifs Discovery"
    """

    def __init__(self, network: DiGraph, isomorphic_mapping: dict):
        super().__init__(network, isomorphic_mapping)
        self.undirected_network = nx.Graph(network)
        self.unique = set()  # unique sub graphs visited

    def __is_unique(self, sub_graph: DiGraph) -> bool:
        return graph_to_hashed_graph(sub_graph) not in self.unique

    def __extend_sub_graphs(self, sub_graph: set, extension: set, v: int):
        if len(sub_graph) == self.k:
            graph = nx.induced_subgraph(self.network, list(sub_graph))
            if self.__is_unique(graph):
                self.unique.add(graph_to_hashed_graph(graph))
                self._inc_count_w_canonical_label(graph)
        else:
            while len(extension) > 0:
                w = random.sample(extension, 1)[0]
                extension.remove(w)

                w_neighbors = set(self.undirected_network.neighbors(w))
                excl_neighbors = w_neighbors.difference(sub_graph)
                v_ext_new = set([u for u in excl_neighbors if u > v])

                new_extension = extension.union(v_ext_new)
                self.__extend_sub_graphs(sub_graph.union({w}), new_extension, v)

    def search_sub_graphs(self, k: int, allow_self_loops: bool) -> SubGraphSearchResult:
        self.fsl = defaultdict(int)
        self.fsl_fully_mapped = defaultdict(list)
        self.k = k
        self.allow_self_loops = allow_self_loops
        self.unique = set()

        for v in list(self.network.nodes):
            self.logger.debug(f'Node: ({v}):')
            v_neighbors = list(self.undirected_network.neighbors(v))
            v_ext = set([u for u in v_neighbors if u > v])
            self.__extend_sub_graphs({v}, v_ext, v)

        self.fsl = dict(sorted(self.fsl.items()))
        return SubGraphSearchResult(fsl=self.fsl, fsl_fully_mapped=self.fsl_fully_mapped)
