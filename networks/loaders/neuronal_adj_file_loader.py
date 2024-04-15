from networks.loaders.network_loader_strategy import NetworkLoaderStrategy
from networks.network import Network


class NeuronalAdjFileLoader(NetworkLoaderStrategy):
    def __init__(self):
        super().__init__()

    def load(self, *args) -> Network:
        """
        adj_file_path: simple txt format: (v1, v2, w) per line
            whereas v1 -> v2, and w is the number of synapses
        neurons_file_path: contains the names of c.elegans neurons
            whereas each line is a name and the index match the adjacency matrix file
        """
        adj_file_path, neurons_file_path = args

        with open(adj_file_path, "r") as f:
            for line in f.readlines():
                v1, v2, w = tuple(line.strip().split())
                self._load_synapse(v1, v2, w)

        with open(neurons_file_path, "r") as f:
            self.neuron_names = [line.strip() for line in f.readlines()]

        return self._copy_network_params()
