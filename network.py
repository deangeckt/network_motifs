import networkx as nx
import matplotlib.pyplot as plt

from utils.config import Config
from utils.simple_logger import Logger


class Network:
    def __init__(self):
        config = Config()
        self.logger = Logger()
        self.graph = nx.DiGraph()

        self.neuron_names = []
        self.amount_of_synapses = 0
        self.synapse_amount_threshold = int(config.get_property('neuronal', 'synapse_amount_threshold'))
        self.participating_neurons = set()

    def load_adj_file(self, file_path: str, is_synapse=False):
        """
        simple txt format: (v1, v2, w) per line
        whereas v1 -> v2, and w either the weight or number of synapses
        """
        with open(file_path, "r") as f:
            for line in f.readlines():
                v1, v2, w = tuple(line.strip().split())
                if is_synapse:
                    if int(w) >= self.synapse_amount_threshold:
                        self.amount_of_synapses += int(w)
                        self.graph.add_edge(int(v1), int(v2))
                        self.participating_neurons.add(int(v1))
                        self.participating_neurons.add(int(v2))
                else:
                    self.graph.add_edge(int(v1), int(v2))

    def load_neurons_file(self, file_path):
        """
        neuronal names of c.elegans
        whereas each line is a name and the index match the adjacency matrix
        :return:
        """
        with open(file_path, "r") as f:
            self.neuron_names = [line.strip() for line in f.readlines()]

    def log_properties(self):
        self.logger.info(f'Network properties:')
        if self.neuron_names:
            self.logger.info(f'  - Neurons: {len(self.neuron_names)}')
            self.logger.info(f'  - Synapses: {self.amount_of_synapses}')
            self.logger.info(f'  - Participating neurons have at least: {self.synapse_amount_threshold} synapses')

        self.logger.info(f'  - Nodes: {len(self.graph)}')
        self.logger.info(f'  - Edges: {len(self.graph.edges)}')
        self.logger.info('')

    def plot(self):
        # TODO: need better plotting tools
        nx.draw_networkx(self.graph, with_labels=True, node_size=600, node_color='lightgreen')
        plt.show()
