from network import Network
from simple_logger import Logger, LogLvl
from itertools import combinations


class Circuits:
    def __init__(self, network: Network, logger: Logger):
        self.network = network
        self.logger = logger

    def count_circuits(self):
        self.__count_self_loops()
        self.__count_mutual_regulation()
        self.__count_fan_outs()
        self.__count_cascades()

        self.logger.info(f"self loops (x -> x): {self.self_loops}")
        self.logger.info(f"mutual regulation (x -> y, y -> x): {self.mutual_regulation}")
        self.logger.info(f"fan outs (x -> y, x -> z): {self.fan_outs}")
        self.logger.info(f"cascades (x -> y, y -> z): {self.cascades}")

    def __count_self_loops(self):
        """
        Counts the number of self loops (x -> x) in the given network.
        O(n)
        """
        self.logger.debug('--- self loops (x -> x) debugging: --- ')

        count = 0
        for node in range(len(self.network.adj_matrix)):
            if self.network.adj_matrix[node, node]:
                self.logger.debug(f'{node} -> {node}')
                count += 1

        self.self_loops = count

    def __count_mutual_regulation(self):
        """
        Counts the number of mutual regulation (x -> y, y -> x)
        O(n^2). (with list representation: o(n*e))
        """
        self.logger.debug('--- mutual regulation (x -> y, y -> x) debugging: --- ')

        count = 0
        for x in range(len(self.network.adj_matrix)):
            for y in range(len(self.network.adj_matrix)):
                if x == y:
                    continue
                if self.network.adj_matrix[x, y] and self.network.adj_matrix[y, x]:
                    self.logger.debug(f'{self.network.mapped_nodes_reverse[x]} '
                                      f'<-> {self.network.mapped_nodes_reverse[y]}')
                    count += 1

        self.mutual_regulation = count

    def __count_cascades(self):
        """
        Counts the number of cascades (x -> y, y -> z) in the given network.
        O(n * out degree^2) using list representation. (with adj matrix: o(n^3))
        """
        self.logger.debug('--- cascades (x -> y, y -> z) debugging: --- ')

        count = 0
        for x, x_neighbors in enumerate(self.network.nodes):
            for y in x_neighbors:
                if x == y:
                    continue
                y_neighbors = self.network.nodes[y]
                for z in y_neighbors:
                    if z == y or z == x:
                        continue
                    self.logger.debug(f'{self.network.mapped_nodes_reverse[x]} -> '
                                      f'{self.network.mapped_nodes_reverse[y]} -> '
                                      f'{self.network.mapped_nodes_reverse[z]}')

                    count += 1

        self.cascades = count

    def __count_fan_outs(self):
        """
        Counts the number of cascades (x -> y, x -> z) in the given network.
        O(n)
        """
        self.logger.debug('--- fan outs (x -> y, x -> z) debugging: --- ')

        count = 0
        for node, neighbors in enumerate(self.network.nodes):
            without_self_neighbors = [n for n in neighbors if n != node]
            n = len(without_self_neighbors)
            if n < 2:
                continue
            count += (n * (n - 1)) / 2

            if self.logger.lvl == LogLvl.debug:
                comb = list(combinations(without_self_neighbors, 2))
                for y_, z_ in comb:
                    x = self.network.mapped_nodes_reverse[node]
                    y = self.network.mapped_nodes_reverse[y_]
                    z = self.network.mapped_nodes_reverse[z_]
                    self.logger.debug(f'{x} -> {y}, {x} -> {z}')

        self.fan_outs = int(count)

        def __count_feed_forward():
            pass
