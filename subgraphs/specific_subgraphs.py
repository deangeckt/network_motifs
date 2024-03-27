from typing import Optional

from networkx import DiGraph

from subgraphs.sub_graphs import SubGraphs
from subgraphs.sub_graphs_utils import MotifName, UniqueSubGraph
from utils.simple_logger import LogLvl
from itertools import combinations


class SpecificSubGraphs(SubGraphs):
    """
    None induced
    """
    def __init__(self, network: DiGraph, search: Optional[list[MotifName]] = None):
        super().__init__(network)
        self.implemented_sub_graphs_search = {MotifName.self_loops: self.__count_self_loops,
                                              MotifName.mutual_regulation: self.__count_mutual_regulation,
                                              MotifName.fan_outs: self.__count_fan_outs,
                                              MotifName.cascades: self.__count_cascades,
                                              MotifName.feed_forwards: self.__count_feed_forward,
                                              MotifName.bi_fan: self.__count_bi_fan
                                              }
        self.search: list[MotifName] = self.implemented_sub_graphs_search.keys() if search is None else search

    def search_sub_graphs(self, k: int) -> dict:
        self.logger.info('Specific Sub Graphs search:')
        res = {}
        for sub_graph in self.search:
            if sub_graph not in self.implemented_sub_graphs_search:
                self.logger.info(f'sub graph {sub_graph} not implemented yet')
                continue
            res[sub_graph] = self.implemented_sub_graphs_search[sub_graph]()

        return res

    def __count_self_loops(self):
        """
        Counts the number of self loops (x -> x) in the given network.
        O(n)
        """
        self.logger.debug('--- self loops (x -> x) debugging: --- ')

        count = 0
        for node in list(self.network.nodes):
            if self.network.has_edge(node, node):
                self.logger.debug(f'{node} -> {node}')
                count += 1

        self.logger.info(f"self loops (x -> x): {count}")
        return count

    def __count_mutual_regulation(self):
        """
        Counts the number of mutual regulation (x -> y, y -> x)
        O(n^2). (with list representation: o(n*e))
        """
        self.logger.debug('--- mutual regulation (x -> y, y -> x) debugging: --- ')

        count = 0
        for x in list(self.network.nodes):
            for y in list(self.network.nodes):
                if x == y:
                    continue
                if self.network.has_edge(x, y) and self.network.has_edge(y, x):
                    self.logger.debug(f'{x} <-> {y}')
                    count += 1

        self.logger.info(f"mutual regulation (x -> y, y -> x): {count}")
        return count

    def __count_cascades(self):
        """
        Counts the number of cascades (x -> y, y -> z) in the given network.
        O(n_i * degree_i^2) using list representation. (with adj matrix: o(n^3))
        """
        self.logger.debug('--- cascades (x -> y, y -> z) debugging: --- ')

        count = 0
        self.debug_hash_ = set()

        for x in list(self.network.nodes):
            x_neighbors = list(self.network.adj[x])
            for y in x_neighbors:
                if x == y:
                    continue
                y_neighbors = list(self.network.adj[y])
                for z in y_neighbors:
                    if z == y or z == x:
                        continue

                    sub_graph = ((x, y), (y, z))
                    sub_graph = UniqueSubGraph(sub_graph)
                    if sub_graph in self.debug_hash_:
                        continue
                    self.debug_hash_.add(sub_graph)
                    self.logger.debug(f'{x} -> {y} -> {z}')
                    count += 1

        self.logger.info(f"cascades (x -> y, y -> z): {count}")
        return count

    def __count_fan_outs(self):
        """
        Counts the number of fan outs (x -> y, x -> z) in the given network.
        O(n)
        """
        self.logger.debug('--- fan outs (x -> y, x -> z) debugging: --- ')

        count = 0
        for x in list(self.network.nodes):
            x_neighbors = list(self.network.adj[x])
            without_self_neighbors = [n for n in x_neighbors if n != x]
            n = len(without_self_neighbors)
            if n < 2:
                continue
            count += (n * (n - 1)) / 2

            if self.logger.lvl == LogLvl.debug:
                comb = list(combinations(without_self_neighbors, 2))
                for y, z in comb:
                    self.logger.debug(f'{x} -> {y}, {x} -> {z}')

        fan_outs = int(count)
        self.logger.info(f"fan outs (x -> y, x -> z): {fan_outs}")
        return fan_outs

    def __count_feed_forward(self):
        """
        Counts the number of feed forward (x -> y, x -> z, y -> z) in the given network.
        O(n * degree^3)
        """
        self.logger.debug('--- feed forwards (x -> y, x -> z, y -> z) debugging: --- ')

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
                    count += 1

        self.logger.info(f"feed forwards (x -> y, x -> z, y -> z): {count}")
        return count

    def __count_bi_fan(self):
        """
        Counts the number of bi fan (k=4) (x -> w, x -> z, y -> w, y -> z) in the given network.
        O(n^2 * degree^2)
        """
        self.logger.debug('--- bi fan (x -> w, x -> z, y -> w, y -> z) debugging: --- ')

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
                        sub_graph = UniqueSubGraph(sub_graph)
                        if sub_graph in hash_:
                            continue
                        hash_.add(sub_graph)
                        self.logger.debug(f'{x} -> {w}, {x} -> {z}, {y} -> {w}, {y} -> {z}')
                        count += 1

        self.logger.info(f"bi fan (x -> w, x -> z, y -> w, y -> z): {count}")
        return count

