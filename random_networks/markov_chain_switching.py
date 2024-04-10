import random

from networkx import DiGraph
from tqdm import tqdm

from random_networks.network_randomizer_abc import NetworkRandomizer
from utils.config import Config


class MarkovChainSwitching(NetworkRandomizer):
    def __init__(self, network: DiGraph):
        super().__init__(network)

        config = Config()

        self.switch_factor = int(config.get_property('random', 'markov_chain_switch_factor'))
        self.markov_chain_num_iterations = network.number_of_edges() * self.switch_factor
        self.success_switch = 0

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

    def __markov_chain(self) -> DiGraph:
        """
        - R. Kannan, P. Tetali, S. Vempala, Random Struct. Algorithms 14, 293 (1999).
        - R. Milo,1, 2 N. Kashtan,2, 3 S. Itzkovitz,1, 2 M. E. J. Newman,4 and U. Alon1,
            "On the uniform generation of random graphs with prescribed degree sequences"
            * Degree constrain is saved
            * Mutual / Double edges (number) is NOT saved
        """
        network: DiGraph = self.network.copy()

        for _ in range(self.markov_chain_num_iterations):
            e1, e2 = random.sample(network.edges, 2)
            x1, y1 = e1
            x2, y2 = e2

            if x1 == x2 or x1 == y2 or x2 == y1 or y1 == y2:
                continue

            if network.has_edge(x1, y2) or network.has_edge(x2, y1):
                continue

            self.success_switch += 1
            self.__switch(network, x1, y1, x2, y2)

        return network
