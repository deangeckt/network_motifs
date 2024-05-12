from abc import ABCMeta, abstractmethod

from networks.network import Network
from utils.simple_logger import Logger
from utils.types import NetworkLoaderArgs


class NetworkLoaderStrategy(metaclass=ABCMeta):
    def __init__(self, args: NetworkLoaderArgs):
        self.logger = Logger()
        self.network = Network(args.synapse_threshold)
        self.args = args

    def _load_synapse(self, v1, v2, num_of_synapse, polarity=None):
        self.network.participating_neurons.add(int(v1))
        self.network.participating_neurons.add(int(v2))
        self.network.amount_of_synapses_in_total += int(num_of_synapse)

        if int(num_of_synapse) >= self.args.synapse_threshold:
            self.network.amount_of_synapses_in_graph += int(num_of_synapse)
            if polarity is not None:
                self.network.graph.add_edge(int(v1), int(v2), polarity=polarity)
            else:
                self.network.graph.add_edge(int(v1), int(v2))

    @abstractmethod
    def load(self, *args):
        pass
