import os
from typing import Optional

import networkx as nx

from networks.loaders.durbin_file_loader import DurbinFileLoader
from networks.loaders.neuronal_polarity_loader import NeuronalPolarityLoader
from networks.loaders.simple_adj_file_loader import SimpleAdjFileLoader
from networks.loaders.worm_wiring_loader import WormWiringLoader
from networks.network import Network
from utils.simple_logger import Logger
from utils.types import NetworkInputType


class NetworkLoader:
    def __init__(self):
        self.logger = Logger()
        self.network = Network()

    def load_network_file(self,
                          input_type: NetworkInputType,
                          file_path: str,
                          sheet_name: Optional[str] = None) -> Network:
        if input_type == NetworkInputType.simple_adj_txt:
            loader = SimpleAdjFileLoader()
        elif input_type == NetworkInputType.worm_wiring_xlsx:
            loader = WormWiringLoader()
        elif input_type == NetworkInputType.polarity_xlsx:
            loader = NeuronalPolarityLoader()
        elif input_type == NetworkInputType.durbin_txt:
            loader = DurbinFileLoader()
        else:
            err = f'Error reading network input file: {input_type}'
            self.logger.info(err)
            raise Exception(err)

        name = os.path.basename(file_path)
        self.logger.info(f'Network name: {name}')

        network = loader.load(file_path, sheet_name)
        network.properties()
        self.network = network
        return network

    def load_graph(self, graph: nx.DiGraph) -> Network:
        network = Network()
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
