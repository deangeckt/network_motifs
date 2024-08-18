from networks.loaders.network_loader_strategy import NetworkLoaderStrategy
from utils.export_import import import_network
from utils.types import NetworkLoaderArgs, NetworkBinaryFile


class BinaryNetworkLoader(NetworkLoaderStrategy):
    def __init__(self, args: NetworkLoaderArgs):
        """
        Loads a custom intersected network
        """
        super().__init__(args)

    def load(self, *args):
        file_path = args[0]
        network_binary: NetworkBinaryFile = import_network(file_path)
        self.network.neuron_names = network_binary['neuron_names']
        self.network.participating_neurons = network_binary['participating_neurons']
        self.network.graph = network_binary['graph']
