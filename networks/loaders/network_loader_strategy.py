from abc import ABCMeta, abstractmethod

from networks.network import Network
from utils.simple_logger import Logger
from utils.types import NetworkLoaderArgs


class NetworkLoaderStrategy(metaclass=ABCMeta):
    def __init__(self, args: NetworkLoaderArgs):
        self.logger = Logger()
        self.network = Network(args.synapse_threshold)
        self.args = args

    def _append_edge(self, v1: int, v2: int, synapse: int, gap: int, polarity=None):
        self.network.participating_neurons.add(v1)
        self.network.participating_neurons.add(v2)

        self.network.amount_of_synapses_in_total += synapse
        self.network.amount_of_gaps_in_total += gap

        if synapse >= self.args.synapse_threshold or gap >= self.args.synapse_threshold:
            self.network.amount_of_synapses_in_graph += synapse
            self.network.amount_of_gaps_in_graph += gap

            if self.network.graph.has_edge(v1, v2):
                edge = self.network.graph.get_edge_data(v1, v2)
                self.network.graph.add_edge(v1,
                                            v2,
                                            synapse=edge['synapse'] + synapse,
                                            gap=edge['gap'] + gap,
                                            polarity=edge['polarity'])
            else:
                self.network.graph.add_edge(v1, v2, synapse=synapse, gap=gap, polarity=polarity)

    @abstractmethod
    def load(self, *args):
        pass
