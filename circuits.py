from network import Network
from utils import Verbose
from itertools import combinations


class Circuits:
    def __init__(self, network: Network, verbose=Verbose.info):
        self.network = network
        self.verbose = verbose

    def count_circuits(self):
        self.__count_self_loops()
        self.__count_mutual_regulation()
        self.__count_fan_outs()
        self.__count_cascades()

        if self.verbose >= Verbose.info:
            print(f"self loops (x -> x): {self.self_loops}")
            print(f"mutual regulation (x -> y, y -> x): {self.mutual_regulation}")
            print(f"fan outs (x -> y, x -> z): {self.fan_outs}")
            print(f"cascades (x -> y, y -> z): {self.cascades}")

    def __count_self_loops(self):
        """
        Counts the number of self loops (x -> x) in the given network.
        O(n)
        """
        if self.verbose >= Verbose.debug:
            print('--- self loops (x -> x) debugging: --- ')

        count = 0
        for node in range(len(self.network.adj_matrix)):
            if self.network.adj_matrix[node, node]:
                if self.verbose >= Verbose.debug:
                    print(f'{node} -> {node}')
                count += 1

        self.self_loops = count

    def __count_mutual_regulation(self):
        """
        Counts the number of mutual regulation (x -> y, y -> x)
        O(n^2). (with list representation: o(n*e))
        """
        if self.verbose >= Verbose.debug:
            print('--- mutual regulation (x -> y, y -> x) debugging: --- ')

        count = 0
        for x in range(len(self.network.adj_matrix)):
            for y in range(len(self.network.adj_matrix)):
                if x == y:
                    continue
                if self.network.adj_matrix[x, y] and self.network.adj_matrix[y, x]:
                    if self.verbose >= Verbose.debug:
                        print(self.network.mapped_nodes_reverse[x], '<->',
                              self.network.mapped_nodes_reverse[y])
                    count += 1

        self.mutual_regulation = count

    def __count_cascades(self):
        """
        Counts the number of cascades (x -> y, y -> z) in the given network.
        O(n * out degree^2) using list representation. (with adj matrix: o(n^3))
        """
        if self.verbose >= Verbose.debug:
            print('--- cascades (x -> y, y -> z) debugging: --- ')

        count = 0
        for x, x_neighbors in enumerate(self.network.nodes):
            for y in x_neighbors:
                if x == y:
                    continue
                y_neighbors = self.network.nodes[y]
                for z in y_neighbors:
                    if z == y or z == x:
                        continue
                    if self.verbose >= Verbose.debug:
                        print(self.network.mapped_nodes_reverse[x], '->',
                              self.network.mapped_nodes_reverse[y], '->',
                              self.network.mapped_nodes_reverse[z])

                    count += 1

        self.cascades = count

    def __count_fan_outs(self):
        """
        Counts the number of cascades (x -> y, x -> z) in the given network.
        O(n)
        """
        if self.verbose >= Verbose.debug:
            print('--- fan outs (x -> y, x -> z) debugging: --- ')

        count = 0
        for node, neighbors in enumerate(self.network.nodes):
            without_self_neighbors = [n for n in neighbors if n != node]
            n = len(without_self_neighbors)
            if n < 2:
                continue
            count += (n * (n - 1)) / 2

            if self.verbose >= Verbose.debug:
                comb = list(combinations(without_self_neighbors, 2))
                for y_, z_ in comb:
                    x = self.network.mapped_nodes_reverse[node]
                    y = self.network.mapped_nodes_reverse[y_]
                    z = self.network.mapped_nodes_reverse[z_]
                    print(f'{x} -> {y}, {x} -> {z}')

        self.fan_outs = int(count)
