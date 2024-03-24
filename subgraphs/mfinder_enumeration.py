from networkx import DiGraph

from subgraphs.sub_graphs import SubGraphs
import networkx as nx

from subgraphs.sub_graphs_utils import get_id, get_sub_graph_from_id
from utils.simple_logger import LogLvl, Logger
from collections import defaultdict


class UniqueSubGraph:
    def __init__(self, sub_graph: tuple):
        self.sub_graph = sub_graph

    def __eq__(self, other):
        return sorted(self.sub_graph) == sorted(other.sub_graph)

    def __hash__(self):
        return hash(tuple(sorted(self.sub_graph)))


class MFinder(SubGraphs):
    """
    MFinder enumeration algorithm
    """

    def __init__(self, network: DiGraph):
        super().__init__(network)
        self.fsl = defaultdict(int)  # frequent sub graph list
        self.fsl_ids = {}
        self.k = -1  # motif size
        self.unique = set()  # unique sub graphs visited
        self.hash_ = set()  # hash for trimming during the backtracking

    def __is_unique(self, sub_graph: tuple) -> bool:
        return UniqueSubGraph(sub_graph) not in self.unique

    def __inc_count_w_canonical_label(self, sub_graph: tuple):
        sub_id = get_id(sub_graph)
        self.logger.debug(f'inc count to motif id: {sub_id}')
        self.logger.debug(sub_graph)

        graph_hash = nx.weisfeiler_lehman_graph_hash(nx.DiGraph(list(sub_graph)), iterations=self.k)
        self.fsl[graph_hash] += 1
        self.fsl_ids[graph_hash] = sub_id

    def __find_sub_graphs_new_edge(self, sub_graph: tuple, edge: tuple):
        if edge in sub_graph:
            return
        new_sub_graph = (*sub_graph, edge)
        if UniqueSubGraph(new_sub_graph) in self.hash_:
            return
        self.__find_sub_graphs(new_sub_graph)

    def __find_sub_graphs(self, sub_graph: tuple):
        graph = nx.DiGraph(list(sub_graph))
        if len(graph) > self.k:
            return
        if len(graph) == self.k and self.__is_unique(sub_graph):
            self.unique.add(UniqueSubGraph(sub_graph))
            self.__inc_count_w_canonical_label(sub_graph)
            if self.k > 2:
                return

        self.hash_.add(UniqueSubGraph(sub_graph))
        for i in list(graph.nodes):
            for k in list(self.network.adj[i]):
                self.__find_sub_graphs_new_edge(sub_graph, (i, k))

            in_edges = [u for u, v in list(self.network.edges) if v == i]  # O(E) - can use adj mat's col for O(N)
            for k in in_edges:
                self.__find_sub_graphs_new_edge(sub_graph, (k, i))

    def search_sub_graphs(self, k: int) -> dict:
        self.fsl = defaultdict(int)
        self.fsl_ids = {}
        self.k = k
        self.unique = set()
        self.hash_ = set()

        self.logger.info('Sub Graphs search:')
        for i, j in list(self.network.edges):
            self.logger.debug(f'({i}, {j})')
            self.__find_sub_graphs(((i, j),))

        fsl_mapped = {self.fsl_ids[hash_]: self.fsl[hash_] for hash_ in self.fsl}
        return fsl_mapped
