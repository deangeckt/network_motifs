import random

from networkx import DiGraph
from tqdm import tqdm

from random_networks.network_randomizer_abc import NetworkRandomizer
from utils.config import Config


class MarkovChainSwitching(NetworkRandomizer):
    def __init__(self, network: DiGraph):
        super().__init__(network)

        config = Config()

        self.q = int(config.get_property('random', 'markov_chain_q'))
        self.markov_chain_num_iterations = network.number_of_edges() * self.q
        self.use_self_loops = config.get_boolean_property('run_args', 'allow_self_loops')

    def generate(self, amount: int) -> list[DiGraph]:
        self.logger.info('\nRandomizer: markov chain switching algo')
        self.logger.info(f'Randomizer: generating {amount} random networks')
        self.logger.info(f'Markov chain q: {self.q}')
        self.logger.info(f'Markov chain iterations: {self.markov_chain_num_iterations}')
        return [self.__markov_chain() for _ in tqdm(range(amount))]

    @staticmethod
    def __switch(network: DiGraph, x1: int, y1: int, x2: int, y2: int):
        network.remove_edge(x1, y1)
        network.remove_edge(x2, y2)
        network.add_edge(x1, y2)
        network.add_edge(x2, y1)

    @staticmethod
    def __can_switch(network: DiGraph, x1: int, y1: int, x2: int, y2: int) -> bool:
        # already exist
        if network.has_edge(x1, y2) or network.has_edge(x2, y1):
            return False

        # double edge
        if network.has_edge(y2, x1) or network.has_edge(y1, x2):
            return False

        return True

    def __markov_chain(self) -> DiGraph:
        """
        - R. Kannan, P. Tetali, S. Vempala, Random Struct. Algorithms 14, 293 (1999).
        - R. Milo,1, 2 N. Kashtan,2, 3 S. Itzkovitz,1, 2 M. E. J. Newman,4 and U. Alon1,
            "On the uniform generation of random graphs with prescribed degree sequences"
            * Degree constrain is saved
            * Mutual edges (number) is saved
        """
        network: DiGraph = self.network.copy()

        for i in range(self.markov_chain_num_iterations):
            x1, x2 = random.sample(list(network.nodes), 2)
            if not network.adj[x1] or not network.adj[x2]:
                continue
            y1 = random.sample(list(network.adj[x1]), 1)[0]
            y2 = random.sample(list(network.adj[x2]), 1)[0]

            if x1 == y2 or x2 == y1:
                continue

            # handle mutual edge (switch only with a different double edge)
            if network.has_edge(y1, x1) and network.has_edge(y2, x2):
                if self.__can_switch(network, x1, y1, x2, y2) and self.__can_switch(network, y1, x1, y2, x2):
                    self.__switch(network, x1, y1, x2, y2)
                    self.__switch(network, y1, x1, y2, x2)
            # handle single edge
            else:
                if self.__can_switch(network, x1, y1, x2, y2):
                    self.__switch(network, x1, y1, x2, y2)

        return network
