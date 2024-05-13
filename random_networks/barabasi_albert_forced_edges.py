import math
import networkx as nx
from networkx import DiGraph
from tqdm import tqdm

from networks.network import Network
from random_networks.network_randomizer_abc import NetworkRandomizer


class BarabasiAlbertForcedEdges(NetworkRandomizer):
    """
    Barabási–Albert model with forced (average) number of edges:
    - the generated graph by the model is undirected one with |E| = n * m, thus m = ciel(e|/n)
    and remove the extra edges
    - We will convert it to directed by choosing randomly a direction for each edge
        * Degree constrain is NOT saved
        * Mutual / Double edges (number) is NOT saved
        * polarity ratio is saved
    - Note: it's best to keep the actual values of the network nodes as a sequence: (0 .. n-1)
            i.e.: a network with only "4->5" edge might generate something like "0 -> 1" - those different names
            might later be counted differently by the sub graph counting algorithm.
    """

    def __init__(self, network: Network):
        super().__init__(network)

        graph = network.graph
        n = len(graph)
        e = len(graph.edges)
        m = math.ceil(e / n)

        self.n = n
        self.e = e
        self.m = m

        self.generate_foo = self.__generate_with_polarity if network.use_polarity else self.__generate

    def generate(self, amount: int) -> list[DiGraph]:
        self.logger.info(f'm (# of edges to attach): {self.m}')

        random_networks = [self.generate_foo() for _ in tqdm(range(amount))]
        self._log_avg_num_of_generated_edges(random_networks, amount)

        return random_networks

    def __generate(self) -> DiGraph:
        undirected_graph = nx.barabasi_albert_graph(n=self.n, m=self.m)
        remove_edges_amount = len(undirected_graph.edges) - self.e
        self._remove_random_edges(undirected_graph, remove_edges_amount)

        # this kills the clustering coefficients properties of the original network
        return self._assign_direction(undirected_graph)

    def __generate_with_polarity(self) -> DiGraph:
        graph = self.__generate()
        self._assign_polarity(graph)
        return graph
