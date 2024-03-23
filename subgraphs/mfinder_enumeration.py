from networkx import DiGraph

from subgraphs.sub_graphs import SubGraphs
import networkx as nx

from subgraphs.sub_graphs_utils import get_id, get_sub_graph_from_id
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
        self.hash_ = set()  # hash for trimming during the backtracking

    def __is_unique(self, sub_graph: tuple) -> bool:
        return UniqueSubGraph(sub_graph) not in self.unique

    def __inc_count_w_canonical_label(self, sub_graph: tuple):
        # TODO: improve run time with hashing as a key - not looping on fsl. check correctness
        sub_id = get_id(sub_graph)
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
            other_graph = get_sub_graph_from_id(decimal=id_, k=self.k)
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
    k = 2
    g = nx.DiGraph([(0, 1), (1, 0), (0, 2), (1, 2), (2, 1), (2, 0), (0, 0)])
    g = nx.DiGraph([(0, 1), (1, 0), (0, 0)])
    mfinder = MFinder(g)
    sub_graphs = mfinder.search_sub_graphs(k=k)
    total = 0
    for sub_id in sub_graphs:
        sub_graph = get_sub_graph_from_id(sub_id, k=k)
        amount = sub_graphs[sub_id]
        total += amount
        print()
        print(nx.adjacency_matrix(sub_graph).todense())
        print(f'Amount: {amount}')

    print('\nmotif amount', len(sub_graphs))
    print('total', total)
