import os

import networkx as nx

from networks.loaders.binary_network_loader import BinaryNetworkLoader
from networks.loaders.durbin_file_loader import DurbinFileLoader
from networks.loaders.multilayer_loader import MultilayerConnectomeLoader
from networks.loaders.neuronal_polarity_loader import NeuronalPolarityLoader
from networks.loaders.scipy_sparse_loader import ScipySparseLoader
from networks.loaders.simple_adj_file_loader import SimpleAdjFileLoader
from networks.loaders.worm_wiring_loader import WormWiringLoader
from networks.network import Network
from utils.simple_logger import Logger
from utils.types import NetworkInputType, NetworkLoaderArgs


class NetworkLoader:
    def __init__(self, args: NetworkLoaderArgs):
        self.logger = Logger()
        self.network = Network(synapse_threshold=args.synapse_threshold)
        self.args = args

    def load_network_file(self, input_type: NetworkInputType, file_path: str) -> Network:
        if input_type == NetworkInputType.simple_adj_txt:
            loader = SimpleAdjFileLoader(self.args)
        elif input_type == NetworkInputType.worm_wiring_xlsx:
            loader = WormWiringLoader(self.args)
        elif input_type == NetworkInputType.polarity_xlsx:
            loader = NeuronalPolarityLoader(self.args)
        elif input_type == NetworkInputType.durbin_txt:
            loader = DurbinFileLoader(self.args)
        elif input_type == NetworkInputType.multilayer:
            loader = MultilayerConnectomeLoader(self.args)
        elif input_type == NetworkInputType.binary_network:
            loader = BinaryNetworkLoader(self.args)
        elif input_type == NetworkInputType.scipy_sparse:
            loader = ScipySparseLoader(self.args)
        else:
            err = f'Error reading network input file: {input_type}'
            self.logger.info(err)
            raise Exception(err)

        name = os.path.basename(file_path)
        self.logger.info(f'Network file name: {name}')

        loader.load(file_path)

        self.network = loader.network
        self.network.properties()

        return self.network

    def load_graph(self, graph: nx.DiGraph) -> Network:
        network = Network(synapse_threshold=self.args.synapse_threshold)
        network.graph = graph
        network.calc_polarity_ratio()
        network.properties()

        self.network = network
        return network

    def export_to_simple_format_file(self, output_file_path: str):
        """
        Exports the network to a simple file format (v1, v2, w) per line.
        w isn't the weight but the number of synapses
        """
        with open(output_file_path, "w") as f:
            for i, j in self.network.graph.edges:
                f.write(f'{i} {j} 1\n')
