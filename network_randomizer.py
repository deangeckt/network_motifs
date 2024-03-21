import random

from networkx import DiGraph
from tqdm import tqdm

from utils.config import Config
from utils.simple_logger import Logger


class NetworkRandomizer:
    def __init__(self, network: DiGraph):
        self.real_network = network
        self.logger = Logger()

        config = Config()
        self.markov_chain_num_iterations = int(config.get_property('random', 'markov_chain_num_iterations'))

    def generate(self, amount: int) -> list[DiGraph]:
        self.logger.info(f'Randomizer: generating {amount} networks using '
                         f'markov chain with {self.markov_chain_num_iterations} iterations')
        return [self.markov_chain() for _ in tqdm(range(amount))]

    def markov_chain(self) -> DiGraph:
        """
        S11. R. Kannan, P. Tetali, S. Vempala, Random Struct. Algorithms 14, 293 (1999).
        Degree constrain is saved
        """
        network: DiGraph = self.real_network.copy()

        for i in range(self.markov_chain_num_iterations):
            x1, x2 = random.sample(list(network.nodes), 2)
            if not network.adj[x1] or not network.adj[x2]:
                continue
            y1 = random.sample(list(network.adj[x1]), 1)[0]
            y2 = random.sample(list(network.adj[x2]), 1)[0]

            if network.has_edge(x1, y2) or network.has_edge(x2, y1):
                continue

            network.remove_edge(x1, y1)
            network.remove_edge(x2, y2)
            network.add_edge(x1, y2)
            network.add_edge(x2, y1)

        return network
