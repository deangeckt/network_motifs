from abc import ABCMeta, abstractmethod

import networkx as nx

from networks.network import Network
from utils.config import Config
from utils.simple_logger import Logger


class NetworkLoaderStrategy(metaclass=ABCMeta):
    def __init__(self):
        config = Config()
        self.logger = Logger()
        self.graph = nx.DiGraph()

        # neurons configuration
        self.neuron_names: list[str] = []
        self.amount_of_synapses_in_graph = 0
        self.amount_of_synapses_in_total = 0
        self.synapse_amount_threshold = int(config.get_property('neuronal', 'synapse_amount_threshold'))
        self.participating_neurons = set()

        # polarity configuration
        self.use_polarity = False

    def _load_synapse(self, v1, v2, num_of_synapse, polarity=None):
        self.participating_neurons.add(int(v1))
        self.participating_neurons.add(int(v2))
        self.amount_of_synapses_in_total += int(num_of_synapse)

        if int(num_of_synapse) >= self.synapse_amount_threshold:
            self.amount_of_synapses_in_graph += int(num_of_synapse)
            if polarity is not None:
                self.graph.add_edge(int(v1), int(v2), polarity=polarity)
            else:
                self.graph.add_edge(int(v1), int(v2))

    def _copy_network_params(self) -> Network:
        network = Network()
        network.graph = self.graph

        # neurons configuration
        network.neuron_names = self.neuron_names
        network.amount_of_synapses_in_graph = self.amount_of_synapses_in_graph
        network.amount_of_synapses_in_total = self.amount_of_synapses_in_total
        network.participating_neurons = self.participating_neurons

        # polarity configuration
        network.use_polarity = self.use_polarity
        return network

    @abstractmethod
    def load(self, *args):
        pass
