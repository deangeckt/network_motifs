from networkx import DiGraph

from subgraphs.sub_graphs import SubGraphs
import networkx as nx

from subgraphs.sub_graphs_utils import get_id, UniqueSubGraph
from collections import defaultdict


class MFinderNoneInduced(SubGraphs):
    """
    MFinder enumeration algorithm. none-induced version of the algo:
    paper:  R. Milo, S. Shen-Orr, S. Itzkovitz, N. Kashtan,D. Chklovskii, and U. Alon, “Network motifs: simple
            building blocks of complex networks.” Science, vol. 298, no. 5594, pp. 824–827, October 2002.
    pseudocode: Ribeiro, Pedro and Silva, Fernando and Kaiser, Marcus: "Strategies for Network Motifs Discovery"
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
        graph = nx.DiGraph(list(sub_graph))
        sub_id = get_id(graph)
        self.logger.debug(f'inc count to motif id: {sub_id}')
        self.logger.debug(sub_graph)

        graph_hash = nx.weisfeiler_lehman_graph_hash(graph, iterations=self.k)
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

            in_edges = [u for u, v in self.network.in_edges(nbunch=i)]
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
            self.logger.debug(f'Edge: ({i}, {j}):')
            self.__find_sub_graphs(((i, j),))

        fsl_mapped = {self.fsl_ids[hash_]: self.fsl[hash_] for hash_ in self.fsl}
        return fsl_mapped