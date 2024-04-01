from functools import cache

from networkx import DiGraph

from subgraphs.sub_graphs_abc import SubGraphsABC
import networkx as nx

from subgraphs.sub_graphs_utils import get_id, graph_to_hashed_graph
from collections import defaultdict

from utils.config import Config


class MFinderInduced(SubGraphsABC):
    """
    MFinder enumeration algorithm.
    paper:  R. Milo, S. Shen-Orr, S. Itzkovitz, N. Kashtan,D. Chklovskii, and U. Alon, “Network motifs: simple
            building blocks of complex networks.” Science, vol. 298, no. 5594, pp. 824–827, October 2002.
    pseudocode: Ribeiro, Pedro and Silva, Fernando and Kaiser, Marcus: "Strategies for Network Motifs Discovery"
    """

    def __init__(self, network: DiGraph, isomorphic_mapping: dict):
        super().__init__(network, isomorphic_mapping)
        config = Config()
        self.map_sub_graphs = config.get_boolean_property('run_args', 'log_all_motif_sub_graphs')

        self.fsl = defaultdict(int)  # frequent sub graph list
        self.k = -1  # motif size
        self.unique = set()  # unique sub graphs visited
        self.hash_ = set()  # hash for trimming during the backtracking
        self.fsl_fully_mapped = defaultdict(list)  # in case we want list of all sub graphs

    def __is_unique(self, sub_graph: DiGraph) -> bool:
        return graph_to_hashed_graph(sub_graph) not in self.unique

    def __inc_count_w_canonical_label(self, sub_graph: DiGraph):
        sub_id = get_id(sub_graph)
        self.logger.debug(f'inc count to motif id: {sub_id}')
        self.logger.debug(list(sub_graph.edges))

        if sub_id not in self.isomorphic_mapping:
            return

        sub_id_isomorphic_representative = self.isomorphic_mapping[sub_id]
        self.fsl[sub_id_isomorphic_representative] += 1

        if self.map_sub_graphs:
            self.fsl_fully_mapped[sub_id_isomorphic_representative].append(list(sub_graph.edges))

    def __find_sub_graphs_new_edge(self, sub_graph: frozenset, k: int):
        if k in sub_graph:
            return
        new_sub_graph = sub_graph.union({k})
        if new_sub_graph in self.hash_:
            return
        self.__find_sub_graphs(new_sub_graph)

    @cache
    def __find_sub_graphs(self, sub_graph: frozenset):
        graph = nx.induced_subgraph(self.network, list(sub_graph))
        if len(graph) > self.k:
            return
        if len(graph) == self.k and self.__is_unique(graph):
            self.unique.add(graph_to_hashed_graph(graph))
            self.__inc_count_w_canonical_label(graph)
            if self.k > 2:
                return

        self.hash_.add(sub_graph)
        for i in list(graph.nodes):
            for k in list(self.network.adj[i]):
                self.__find_sub_graphs_new_edge(sub_graph, k)

            in_edges = [u for u, v in self.network.in_edges(nbunch=i)]
            for k in in_edges:
                self.__find_sub_graphs_new_edge(sub_graph, k)

    def search_sub_graphs(self, k: int) -> dict:
        self.fsl = defaultdict(int)
        self.k = k
        self.unique = set()
        self.hash_ = set()
        self.fsl_fully_mapped = defaultdict(list)

        for i, j in list(self.network.edges):
            self.logger.debug(f'Edge: ({i}, {j}):')
            self.__find_sub_graphs(frozenset({i, j}))

        return dict(sorted(self.fsl.items()))

    def get_sub_graphs_fully_mapped(self) -> dict:
        return self.fsl_fully_mapped
