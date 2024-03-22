import numpy as np
from networkx import DiGraph

from subgraphs.sub_graphs import SubGraphs
import networkx as nx

from utils.simple_logger import LogLvl, Logger


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
        self.fsl = {}  # frequent sub graph list
        self.k = -1  # motif size
        self.unique = set()  # unique sub graphs visited
        self.hash_ = set()

    @staticmethod
    def __get_id(sub_graph: tuple) -> int:
        graph = nx.DiGraph(list(sub_graph))
        adj_mat = nx.adjacency_matrix(graph).todense()
        vec = adj_mat.flatten()
        decimal = 0
        for i, bit in enumerate(vec):
            decimal += bit * (2 ** i)
        return decimal

    def get_sub_graph_from_id(self, decimal: int) -> DiGraph:
        bin_digits = [int(d) for d in str(bin(decimal))[2:]]
        pad_amount = self.k ** 2 - len(bin_digits)
        padding = pad_amount * [0]
        bin_digits = padding + bin_digits
        bin_digits.reverse()
        adj_mat = np.array(bin_digits).reshape(self.k, self.k)
        return nx.DiGraph(adj_mat)

    def __is_unique(self, sub_graph: tuple) -> bool:
        return UniqueSubGraph(sub_graph) not in self.unique

    def __inc_count_w_canonical_label(self, sub_graph: tuple):
        sub_id = self.__get_id(sub_graph)
        self.logger.debug(f'inc count to motif id: {sub_id}')
        self.logger.debug(sub_graph)
        if len(self.fsl) == 0:
            self.fsl[sub_id] = 1
            return
        if sub_id in self.fsl:
            self.fsl[sub_id] += 1
            return

        graph = nx.DiGraph(list(sub_graph))
        found = False
        for id_ in self.fsl:
            other_graph = self.get_sub_graph_from_id(id_)
            if nx.is_isomorphic(graph, other_graph):
                self.fsl[id_] += 1
                found = True
                break

        if not found:
            self.fsl[sub_id] = 1

    def __find_sub_graphs_new_edge(self, sub_graph: tuple, edge: tuple):
        if edge in sub_graph:
            return
        new_sub_graph = (*sub_graph, edge)
        if UniqueSubGraph(new_sub_graph) in self.hash_:
            return
        self.__find_sub_graphs(new_sub_graph)

    def __find_sub_graphs(self, sub_graph: tuple):
        graph = nx.DiGraph(list(sub_graph))
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
        self.fsl = {}
        self.k = k
        self.unique = set()
        self.hash_ = set()

        self.logger.info('Sub Graphs search:')
        for i, j in list(self.network.edges):
            self.__find_sub_graphs(((i, j),))

        return self.fsl


if __name__ == "__main__":
    logger = Logger(LogLvl.debug)

    g = nx.DiGraph([(0, 1), (1, 0), (0, 2), (1, 2), (2, 1), (2, 0), (0, 0)])
    g = nx.DiGraph([(0, 1), (1, 0), (0, 0)])
    mfinder = MFinder(g)
    sub_graphs = mfinder.search_sub_graphs(k=2)
    total = 0
    for sub_id in sub_graphs:
        sub_graph = mfinder.get_sub_graph_from_id(sub_id)
        amount = sub_graphs[sub_id]
        total += amount
        print()
        print(nx.adjacency_matrix(sub_graph).todense())
        print(f'Amount: {amount}')

    print('\nmotif amount', len(sub_graphs))
    print('total', total)
