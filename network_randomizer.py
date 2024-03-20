import random

from tqdm import tqdm

from network import Network
from utils.config import Config
from utils.simple_logger import Logger
import copy


class NetworkRandomizer:

    def __init__(self, network: Network):
        self.real_network = network
        self.logger = Logger()
        config = Config()

        self.markov_chain_num_iterations = int(config.get_property('random', 'markov_chain_num_iterations'))

    def generate(self, amount: int) -> list[Network]:
        self.logger.info(f'randomizer: markov chain with {self.markov_chain_num_iterations} iterations')
        return [self.markov_chain() for _ in tqdm(range(amount))]

    def markov_chain(self) -> Network:
        """
        S11. R. Kannan, P. Tetali, S. Vempala, Random Struct. Algorithms 14, 293 (1999).
        Degree constrain is saved
        """
        network = copy.deepcopy(self.real_network)
        network_size = len(network.nodes)

        for i in range(self.markov_chain_num_iterations):
            x1, x2 = random.sample(list(range(network_size)), 2)
            if not network.nodes[x1] or not network.nodes[x2]:
                continue
            y1 = random.sample(network.nodes[x1], 1)[0]
            y2 = random.sample(network.nodes[x2], 1)[0]

            if network.adj_matrix[x1, y2] or network.adj_matrix[x2, y1]:
                continue

            network.adj_matrix[x1, y1] = 0
            network.adj_matrix[x2, y2] = 0

            network.adj_matrix[x1, y2] = 1
            network.adj_matrix[x2, y1] = 1

            network.nodes[x1].remove(y1)
            network.nodes[x2].remove(y2)

            network.nodes[x1].append(y2)
            network.nodes[x2].append(y1)

        return network
