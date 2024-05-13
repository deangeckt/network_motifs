import networkx as nx
from networkx import DiGraph
from tqdm import tqdm

from networks.network import Network
from random_networks.network_randomizer_abc import NetworkRandomizer


class ErdosRenyiForcedEdges(NetworkRandomizer):
    """
    Erdos Renyi algorithm with forced (average) number of edges:
    - there are # n choose 2 possible edges. (times 2 for directed graphs)
    - regular ER graph has an average of (n choose 2 * p) edges; n(n-1)/2 * p. (times 2 for directed graphs)
    - to generate a graph with average |E| edges (given by the real network), we need to set p = |E| / (n(n-1))
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
        p = e / (n * (n - 1))

        self.n = n
        self.e = e
        self.p = p

        self.generate_foo = self.__generate_with_polarity if network.use_polarity else self.__generate

    def generate(self, amount: int) -> list[DiGraph]:
        self.logger.info(f'probability for edge creation: {round(self.p, 4)}')

        random_networks = [self.generate_foo() for _ in tqdm(range(amount))]
        self._log_avg_num_of_generated_edges(random_networks, amount)

        return random_networks

    def __generate(self) -> DiGraph:
        return nx.erdos_renyi_graph(n=self.n, p=self.p, directed=True)

    def __generate_with_polarity(self) -> DiGraph:
        graph = self.__generate()
        self._assign_polarity(graph)
        return graph

