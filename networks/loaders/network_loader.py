from typing import Optional

import networkx as nx

from networks.loaders.durbin_file_loader import DurbinFileLoader
from networks.loaders.neuronal_adj_file_loader import NeuronalAdjFileLoader
from networks.loaders.neuronal_polarity_loader import NeuronalPolarityLoader
from networks.loaders.simple_adj_file_loader import SimpleAdjFileLoader
from networks.loaders.worm_wiring_loader import WormWiringLoader
from networks.network import Network
from utils.simple_logger import Logger


class NetworkLoader:
    def __init__(self):
        self.logger = Logger()
        self.network = Network()

    def load_network_file(
            self,
            name: str,
            adj_file_path: Optional[str] = None,
            neurons_file_path: Optional[str] = None,
            worm_wiring_xlsx_file_path: Optional[str] = None,
            worm_wiring_sheet_name: Optional[str] = None,
            polarity_xlsx_file_path: Optional[str] = None,
            polarity_sheet_name: Optional[str] = None,
            durbin_file_path: Optional[str] = None
    ) -> Network:

        self.logger.info(f'Network name: {name}')

        if adj_file_path and neurons_file_path:
            loader = NeuronalAdjFileLoader()
            network = loader.load(adj_file_path, neurons_file_path)
        elif adj_file_path:
            loader = SimpleAdjFileLoader()
            network = loader.load(adj_file_path)
        elif worm_wiring_xlsx_file_path and worm_wiring_sheet_name:
            loader = WormWiringLoader()
            network = loader.load(worm_wiring_xlsx_file_path, worm_wiring_sheet_name)
        elif polarity_xlsx_file_path and polarity_sheet_name:
            loader = NeuronalPolarityLoader()
            network = loader.load(polarity_xlsx_file_path, polarity_sheet_name)
        elif durbin_file_path:
            loader = DurbinFileLoader()
            network = loader.load(durbin_file_path)
        else:
            err = 'Error reading input file'
            self.logger.info(err)
            raise Exception(err)

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
