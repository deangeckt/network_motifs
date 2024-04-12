from collections import defaultdict

from networkx import DiGraph

from subgraphs.sub_graphs_abc import SubGraphsABC
from subgraphs.sub_graphs_utils import HashedGraph
from itertools import combinations

from utils.types import SubGraphSearchResult


class SpecificSubGraphs(SubGraphsABC):
    """
    None induced
    """

    def __init__(self, network: DiGraph, isomorphic_mapping: dict):
        super().__init__(network, isomorphic_mapping)
        self.fsl = {}
        self.fsl_fully_mapped = defaultdict(list)

        self.two_sub_graphs_search = {
            6: self.__count_mutual_regulation
        }

        self.three_sub_graphs_search = {
            6: self.__count_fan_outs,
            12: self.__count_cascades,
            38: self.__count_feed_forward,
        }

        self.four_sub_graph_search = {
            204: self.__count_bi_fan
        }

        self.sub_graphs_ids_per_k = {
            2: self.two_sub_graphs_search,
            3: self.three_sub_graphs_search,
            4: self.four_sub_graph_search
        }

    def search_sub_graphs(self, k: int) -> SubGraphSearchResult:
        self.fsl = {}
        self.fsl_fully_mapped = defaultdict(list)

        sub_graph_searches = self.sub_graphs_ids_per_k[k]

        for id_ in sub_graph_searches:
            sub_graph_searches[id_](id_)

        return SubGraphSearchResult(fsl=self.fsl, fsl_fully_mapped=self.fsl_fully_mapped)

    def __count_self_loops(self, id_: int):
        """
        Counts the number of self loops (x -> x) in the given network.
        O(n)
        """
        self.logger.debug('--- self loops (x -> x) debugging: --- ')
        self.fsl_fully_mapped[id_] = []

        count = 0
        for node in list(self.network.nodes):
            if self.network.has_edge(node, node):
                self.logger.debug(f'{node} -> {node}')
                self.fsl_fully_mapped[id_].append(((node, node), ))
                count += 1

        self.fsl[id_] = count

    def __count_mutual_regulation(self, id_: int):
        """
        Counts the number of mutual regulation (x -> y, y -> x)
        O(n^2). (with list representation: o(n*e))
        """
        self.logger.debug('--- mutual regulation (x -> y, y -> x) debugging: --- ')
        self.fsl_fully_mapped[id_] = []

        count = 0
        for x in list(self.network.nodes):
            for y in list(self.network.nodes):
                if x == y:
                    continue
                if self.network.has_edge(x, y) and self.network.has_edge(y, x):
                    self.logger.debug(f'{x} <-> {y}')
                    self.fsl_fully_mapped[id_].append(((x, y),))
                    count += 1

        self.fsl[id_] = count

    def __count_cascades(self, id_: int):
        """
        Counts the number of cascades (x -> y, y -> z) in the given network.
        O(n_i * degree_i^2) using list representation. (with adj matrix: o(n^3))
        """
        self.logger.debug('--- cascades (x -> y, y -> z) debugging: --- ')
        self.fsl_fully_mapped[id_] = []

        count = 0

        for x in list(self.network.nodes):
            x_neighbors = list(self.network.adj[x])
            for y in x_neighbors:
                if x == y:
                    continue
                y_neighbors = list(self.network.adj[y])
                for z in y_neighbors:
                    if z == y or z == x:
                        continue

                    self.logger.debug(f'{x} -> {y} -> {z}')
                    self.fsl_fully_mapped[id_].append(((x, y), (y, z)))
                    count += 1

        self.fsl[id_] = count

    def __count_fan_outs(self, id_: int):
        """
        Counts the number of fan outs (x -> y, x -> z) in the given network.
        O(n)
        """
        self.logger.debug('--- fan outs (x -> y, x -> z) debugging: --- ')
        self.fsl_fully_mapped[id_] = []

        count = 0
        for x in list(self.network.nodes):
            x_neighbors = list(self.network.adj[x])
            without_self_neighbors = [n for n in x_neighbors if n != x]
            n = len(without_self_neighbors)
            if n < 2:
                continue
            count += (n * (n - 1)) / 2

            comb = list(combinations(without_self_neighbors, 2))
            for y, z in comb:
                self.logger.debug(f'{x} -> {y}, {x} -> {z}')
                self.fsl_fully_mapped[id_].append(((x, y), (x, z)))

        self.fsl[id_] = int(count)

    def __count_feed_forward(self, id_: int):
        """
        Counts the number of feed forward (x -> y, x -> z, y -> z) in the given network.
        O(n * degree^3)
        """
        self.logger.debug('--- feed forwards (x -> y, x -> z, y -> z) debugging: --- ')
        self.fsl_fully_mapped[id_] = []

        count = 0
        for x in list(self.network.nodes):
            x_neighbors = list(self.network.adj[x])
            for y in x_neighbors:
                if x == y:
                    continue
                y_neighbors = list(self.network.adj[y])
                for z in y_neighbors:
                    if z == y or z == x:
                        continue
                    if z not in x_neighbors:
                        continue
                    self.logger.debug(f'{x} -> {y}, {x} -> {z}, {y} -> {z}')
                    self.fsl_fully_mapped[id_].append(((x, y), (x, z), (y, z)))
                    count += 1

        self.fsl[id_] = count

    def __count_bi_fan(self, id_: int):
        """
        Counts the number of bi fan (k=4) (x -> w, x -> z, y -> w, y -> z) in the given network.
        O(n^2 * degree^2)
        """
        self.logger.debug('--- bi fan (x -> w, x -> z, y -> w, y -> z) debugging: --- ')
        self.fsl_fully_mapped[id_] = []

        count = 0
        nodes = list(self.network.nodes)
        nodes.sort()
        N = len(nodes)
        hash_ = set()
        for i in range(N - 1):
            x = nodes[i]
            x_neighbors = list(self.network.adj[x])
            x_without_self_neighbors = [n for n in x_neighbors if n != x]
            x_without_self_neighbors.sort()
            for j in range(1, N):
                y = nodes[j]
                if x == y:
                    continue
                y_neighbors = list(self.network.adj[y])
                y_without_self_neighbors = [n for n in y_neighbors if n != y]
                y_without_self_neighbors.sort()
                y_comb = list(combinations(y_without_self_neighbors, 2))
                x_comb = list(combinations(x_without_self_neighbors, 2))

                for x_wz in x_comb:
                    if x_wz in y_comb:
                        w, z = x_wz
                        sub_graph = ((x, w), (x, z), (y, w), (y, z))
                        sub_graph = HashedGraph(sub_graph)
                        if sub_graph in hash_:
                            continue
                        hash_.add(sub_graph)
                        self.logger.debug(f'{x} -> {w}, {x} -> {z}, {y} -> {w}, {y} -> {z}')
                        self.fsl_fully_mapped[id_].append(sub_graph)
                        count += 1

        self.fsl[id_] = count
