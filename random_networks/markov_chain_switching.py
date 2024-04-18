import random

from networkx import DiGraph
from tqdm import tqdm

from networks.network import Network
from random_networks.network_randomizer_abc import NetworkRandomizer
from utils.config import Config


class MarkovChainSwitching(NetworkRandomizer):
    """
    - R. Kannan, P. Tetali, S. Vempala, Random Struct. Algorithms 14, 293 (1999).
    - R. Milo,1, 2 N. Kashtan,2, 3 S. Itzkovitz,1, 2 M. E. J. Newman,4 and U. Alon1,
        "On the uniform generation of random graphs with prescribed degree sequences"
        * Degree constrain is saved
        * Mutual / Double edges (number) is NOT saved
        * polarity ratio is saved
    """

    def __init__(self, network: Network):
        super().__init__(network)

        config = Config()

        self.switch_factor = int(config.get_property('random', 'markov_chain_switch_factor'))
        self.markov_chain_num_iterations = network.graph.number_of_edges() * self.switch_factor
        self.success_switch = 0

        self.switch_foo = self.__switch_with_polarity if network.use_polarity else self.__switch

    def generate(self, amount: int) -> list[DiGraph]:
        self.logger.info('\nRandomizer: using markov chain switching algorithm')
        self.logger.info(f'Randomizer: generating {amount} random networks')
        self.logger.info(f'Markov chain switch factor: {self.switch_factor}')
        self.logger.info(f'Markov chain iterations: {self.markov_chain_num_iterations}')
        random_networks = [self.__markov_chain() for _ in tqdm(range(amount))]
        self.logger.info(
            f'Markov chain success ratio: {round(self.success_switch / (self.markov_chain_num_iterations * amount), 5)}')
        return random_networks

    @staticmethod
    def __switch(network: DiGraph, x1: int, y1: int, x2: int, y2: int):
        network.remove_edge(x1, y1)
        network.remove_edge(x2, y2)
        network.add_edge(x1, y2)
        network.add_edge(x2, y1)

    @staticmethod
    def __switch_with_polarity(network: DiGraph, x1: int, y1: int, x2: int, y2: int):
        e1_polarity = network.edges[x1, y1]['polarity']
        e2_polarity = network.edges[x2, y2]['polarity']

        network.remove_edge(x1, y1)
        network.remove_edge(x2, y2)

        network.add_edge(x1, y2, polarity=e1_polarity)
        network.add_edge(x2, y1, polarity=e2_polarity)

    def __markov_chain(self) -> DiGraph:
        graph: DiGraph = self.network.graph.copy()

        for _ in range(self.markov_chain_num_iterations):
            e1, e2 = random.sample(graph.edges, 2)
            x1, y1 = e1
            x2, y2 = e2

            if x1 == x2 or x1 == y2 or x2 == y1 or y1 == y2:
                continue

            if graph.has_edge(x1, y2) or graph.has_edge(x2, y1):
                continue

            self.success_switch += 1
            self.switch_foo(graph, x1, y1, x2, y2)

        return graph
