import random

from networkx import DiGraph

from subgraphs.sub_graphs_abc import SubGraphsABC
import networkx as nx

from subgraphs.sub_graphs_utils import get_id, graph_to_hashed_graph
from collections import defaultdict

from utils.types import SubGraphSearchResult


class FanmodESU(SubGraphsABC):
    """
    FanMOD / ESU
    pseudocode: Ribeiro, Pedro and Silva, Fernando and Kaiser, Marcus: "Strategies for Network Motifs Discovery"
    """

    def __init__(self, network: DiGraph, isomorphic_mapping: dict):
        super().__init__(network, isomorphic_mapping)
        self.undirected_network = nx.Graph(network)
        self.fsl = defaultdict(int)
        self.fsl_fully_mapped = defaultdict(list)

        self.k = -1  # motif size
        self.unique = set()  # unique sub graphs visited

    def __is_unique(self, sub_graph: DiGraph) -> bool:
        return graph_to_hashed_graph(sub_graph) not in self.unique

    def __inc_count_w_canonical_label(self, sub_graph: DiGraph):
        sub_id = get_id(sub_graph)
        if sub_id not in self.isomorphic_mapping:
            return

        self.logger.debug(f'inc count to motif id: {sub_id}')
        self.logger.debug(f'{list(sub_graph.edges)}\n')

        sub_id_isomorphic_representative = self.isomorphic_mapping[sub_id]
        self.fsl[sub_id_isomorphic_representative] += 1
        self.fsl_fully_mapped[sub_id_isomorphic_representative].append(tuple(list(sub_graph.edges)))

    def __extend_sub_graphs(self, sub_graph: set, extension: set, v: int):
        if len(sub_graph) == self.k:
            graph = nx.induced_subgraph(self.network, list(sub_graph))
            if self.__is_unique(graph):
                self.unique.add(graph_to_hashed_graph(graph))
                self.__inc_count_w_canonical_label(graph)
        else:
            extension_set = set(extension)
            while len(extension_set) > 0:
                w = random.sample(extension_set, 1)[0]
                extension_set.remove(w)

                w_neighbors = set(self.undirected_network.neighbors(w))
                excl_neighbors = w_neighbors.difference(sub_graph)
                v_ext_new = set([u for u in excl_neighbors if u > v])

                new_extension = extension_set.union(v_ext_new)
                self.__extend_sub_graphs(sub_graph.union({w}), new_extension, v)

    def search_sub_graphs(self, k: int) -> SubGraphSearchResult:
        self.fsl = defaultdict(int)
        self.fsl_fully_mapped = defaultdict(list)
        self.k = k
        self.unique = set()

        for v in list(self.network.nodes):
            self.logger.debug(f'Node: ({v}):')
            v_neighbors = list(self.undirected_network.neighbors(v))
            v_ext = set([u for u in v_neighbors if u > v])
            self.__extend_sub_graphs({v}, v_ext, v)

        self.fsl = dict(sorted(self.fsl.items()))
        return SubGraphSearchResult(fsl=self.fsl, fsl_fully_mapped=self.fsl_fully_mapped)
