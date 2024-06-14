from functools import cache

from networkx import DiGraph

from subgraphs.sub_graphs_abc import SubGraphsABC
import networkx as nx

from utils.sub_graphs import HashedGraph
from collections import defaultdict

from utils.types import SubGraphSearchResult


class MFinderNoneInduced(SubGraphsABC):
    """
    MFinder enumeration algorithm. none-induced version of the algo:
    paper:  R. Milo, S. Shen-Orr, S. Itzkovitz, N. Kashtan,D. Chklovskii, and U. Alon, “Network motifs: simple
            building blocks of complex networks.” Science, vol. 298, no. 5594, pp. 824–827, October 2002.
    pseudocode: Ribeiro, Pedro and Silva, Fernando and Kaiser, Marcus: "Strategies for Network Motifs Discovery"
    """

    def __init__(self, network: DiGraph, isomorphic_mapping: dict):
        super().__init__(network, isomorphic_mapping)
        self.unique = set()  # unique sub graphs visited
        self.hash_ = set()  # hash for trimming during the backtracking

    def __is_unique(self, sub_graph: tuple) -> bool:
        return HashedGraph(sub_graph) not in self.unique

    def __find_sub_graphs_new_edge(self, sub_graph: tuple, edge: tuple):
        if edge in sub_graph:
            return
        new_sub_graph = (*sub_graph, edge)
        if HashedGraph(new_sub_graph) in self.hash_:
            return
        self.__find_sub_graphs(new_sub_graph)

    @cache
    def __find_sub_graphs(self, sub_graph: tuple):
        graph = nx.DiGraph(list(sub_graph))
        if len(graph) > self.k:
            return
        if len(graph) == self.k and self.__is_unique(sub_graph):
            self.unique.add(HashedGraph(sub_graph))
            self._inc_count_w_canonical_label(nx.DiGraph(list(sub_graph)))
            if self.k > 2:
                return

        self.hash_.add(HashedGraph(sub_graph))
        for i in list(graph.nodes):
            for k in list(self.network.adj[i]):
                self.__find_sub_graphs_new_edge(sub_graph, (i, k))

            in_edges = [u for u, v in self.network.in_edges(nbunch=i)]
            for k in in_edges:
                self.__find_sub_graphs_new_edge(sub_graph, (k, i))

    def search_sub_graphs(self, k: int, allow_self_loops: bool) -> SubGraphSearchResult:
        self.fsl = defaultdict(int)
        self.k = k
        self.allow_self_loops = allow_self_loops
        self.unique = set()
        self.hash_ = set()
        self.fsl_fully_mapped = defaultdict(list)

        for i, j in list(self.network.edges):
            self.logger.debug(f'Edge: ({i}, {j}):')
            self.__find_sub_graphs(((i, j),))

        self.fsl = dict(sorted(self.fsl.items()))
        return SubGraphSearchResult(fsl=self.fsl, fsl_fully_mapped=self.fsl_fully_mapped)
