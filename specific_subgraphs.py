from networkx import DiGraph
from utils.simple_logger import Logger, LogLvl
from itertools import combinations


class SpecificSubGraphs:
    def __init__(self, network: DiGraph, search=None):
        self.network = network
        self.logger = Logger()
        self.implemented_sub_graphs_search = {'self_loops': self.__count_self_loops,
                                              'mutual_regulation': self.__count_mutual_regulation,
                                              'fan_outs': self.__count_fan_outs,
                                              'cascades': self.__count_cascades,
                                              'feed_forwards': self.__count_feed_forward
                                              }
        self.search = self.implemented_sub_graphs_search.keys() if search is None else search

    def count_sub_graphs(self) -> dict:
        self.logger.info('Specific Sub Graphs:')

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
                    count += 1

        self.logger.info(f"fan outs (x -> y, x -> z): {count}")
        return count

    def __count_fan_outs(self):
        """
        Counts the number of cascades (x -> y, x -> z) in the given network.
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
        self.logger.info(f"cascades (x -> y, y -> z): {fan_outs}")
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
