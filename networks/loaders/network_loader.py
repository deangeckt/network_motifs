import os
from typing import Optional, Union

import networkx as nx

from networks.loaders.durbin_file_loader import DurbinFileLoader
from networks.loaders.neuronal_polarity_loader import NeuronalPolarityLoader
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

    def load_network_file(self,
                          input_type: NetworkInputType,
                          file_path: str,
                          sheet_name: Optional[str] = None) -> Network:
        if input_type == NetworkInputType.simple_adj_txt:
            loader = SimpleAdjFileLoader(self.args)
        elif input_type == NetworkInputType.worm_wiring_xlsx:
            loader = WormWiringLoader(self.args)
        elif input_type == NetworkInputType.polarity_xlsx:
            loader = NeuronalPolarityLoader(self.args)
        elif input_type == NetworkInputType.durbin_txt:
            loader = DurbinFileLoader(self.args)
        else:
            err = f'Error reading network input file: {input_type}'
            self.logger.info(err)
            raise Exception(err)

        name = os.path.basename(file_path)
        self.logger.info(f'Network file name: {name}')
        if sheet_name:
            self.logger.info(f'Sheet name: {sheet_name}')

        loader.load(file_path, sheet_name)

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

    def export_to_gephi(self, output_file_path: str):
        mapping = {i: n for i, n in enumerate(self.network.neuron_names)}
        g = nx.relabel_nodes(self.network.graph, mapping) if mapping else self.network.graph
        nx.write_gexf(g, output_file_path)
