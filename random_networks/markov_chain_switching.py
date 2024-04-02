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

    def __markov_chain(self) -> DiGraph:
        """
        S11. R. Kannan, P. Tetali, S. Vempala, Random Struct. Algorithms 14, 293 (1999).
        Degree constrain is saved
        """
        network: DiGraph = self.network.copy()

        for i in range(self.markov_chain_num_iterations):
            x1, x2 = random.sample(list(network.nodes), 2)
            if not network.adj[x1] or not network.adj[x2]:
                continue
            y1 = random.sample(list(network.adj[x1]), 1)[0]
            y2 = random.sample(list(network.adj[x2]), 1)[0]

            if network.has_edge(x1, y2) or network.has_edge(x2, y1):
                continue

            if not self.use_self_loops and (x1 == y2 or x2 == y1):
                continue

            network.remove_edge(x1, y1)
            network.remove_edge(x2, y2)
            network.add_edge(x1, y2)
            network.add_edge(x2, y1)

        return network
