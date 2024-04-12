from typing import Union

import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from networkx import DiGraph

from utils.config import Config
from utils.simple_logger import Logger


class Network:
    def __init__(self):
        config = Config()
        self.logger = Logger()
        self.graph = nx.DiGraph()

        # general configuration
        self.plot = config.get_boolean_property('run_args', 'plot_properties')

        # neurons configuration
        self.neuron_names: list[str] = []
        self.amount_of_synapses_in_graph = 0
        self.amount_of_synapses_in_total = 0
        self.synapse_amount_threshold = int(config.get_property('neuronal', 'synapse_amount_threshold'))
        self.participating_neurons = set()
        self.participating_nodes = set()

        # polarity configuration
        self.p_src_col = int(config.get_property('polarity', 'src_col'))
        self.p_tar_col = int(config.get_property('polarity', 'tar_col'))
        self.p_weight_col = int(config.get_property('polarity', 'weight_col'))
        self.p_edge_type_col = int(config.get_property('polarity', 'edge_type_col'))
        self.p_polarity_col = int(config.get_property('polarity', 'polarity_col'))
        self.p_prim_nt_col = int(config.get_property('polarity', 'prim_nt_col'))

        # polarity options [+, -, no pred, complex]
        self.filter_polarity = config.get_string_list('polarity', 'filter_polarity')
        # primary neurotransmitter options [GABA, Glu, ACh, 0] (0 is an int)
        self.filter_prim_nt = config.get_string_list('polarity', 'filter_prim_nt')

    def load_graph(self, graph: DiGraph):
        self.graph = graph

    def __load_synapse(self, v1, v2, w):
        self.participating_neurons.add(int(v1))
        self.participating_neurons.add(int(v2))
        self.amount_of_synapses_in_total += int(w)

        if int(w) >= self.synapse_amount_threshold:
            self.amount_of_synapses_in_graph += int(w)
            self.graph.add_edge(int(v1), int(v2))
            self.participating_nodes.add(int(v1))
            self.participating_nodes.add(int(v2))

    # TODO: network new folder with loaders (per type of file...) - returning network?
    def load_polarity_neuronal_file(self, xlsx_path: str, sheet_name: str):
        """
        xlsx files from the paper: Fenyves BG, Szilágyi GS, Vassy Z, Sőti C, Csermely P.
        Synaptic polarity and sign-balance prediction using gene expression data in the Caenorhabditis elegans chemical synapse neuronal connectome network
        """
        xls = pd.ExcelFile(xlsx_path)
        df = xls.parse(sheet_name, header=None)

        # filter
        self.logger.info(f'\nFiltering Neurons with polarity: {self.filter_polarity}')
        df = df[df[self.p_polarity_col].isin(self.filter_polarity)]
        self.logger.info(f'Filtering Neurons with primary neurotransmitter: {self.filter_prim_nt}\n')
        df = df[df[self.p_prim_nt_col].isin(self.filter_prim_nt)]

        src_neurons_names = df.iloc[:, self.p_src_col]
        tar_neurons_names = df.iloc[:, self.p_tar_col]
        edge_weights = df.iloc[:, self.p_weight_col]
        # polarity = df.iloc[:, self.p_polarity_col]
        # primary_neurotransmitter = df.iloc[:, self.p_prim_nt_col]

        self.neuron_names = list(set(src_neurons_names) | set(tar_neurons_names))
        neurons_indices = {ss: i for i, ss in enumerate(self.neuron_names)}

        for v1, v2, w in zip(src_neurons_names, tar_neurons_names, edge_weights):
            self.__load_synapse(neurons_indices[v1], neurons_indices[v2], w)

    def load_adj_neuronal_file(self, adj_file_path: str, neurons_file_path: str):
        """
        adj_file_path: simple txt format: (v1, v2, w) per line
            whereas v1 -> v2, and w is the number of synapses
        neurons_file_path: contains the names of c.elegans neurons
            whereas each line is a name and the index match the adjacency matrix file
        """
        with open(adj_file_path, "r") as f:
            for line in f.readlines():
                v1, v2, w = tuple(line.strip().split())
                self.__load_synapse(v1, v2, w)

        with open(neurons_file_path, "r") as f:
            self.neuron_names = [line.strip() for line in f.readlines()]

    def load_adj_file(self, file_path):
        """
        simple txt format: (v1, v2, w) per line
        whereas v1 -> v2, and w is ignored
        """
        with open(file_path, "r") as f:
            for line in f.readlines():
                v1, v2, _ = tuple(line.strip().split())
                self.graph.add_edge(int(v1), int(v2))

    def properties(self):
        self.logger.info(f'Network properties:')
        if self.neuron_names:
            self.logger.info(f'  - Neurons: {len(self.neuron_names)}')
            self.logger.info(f'  - Neurons with a Synapse: {len(self.participating_neurons)}')
            self.logger.info(f'  - Synapses in the network: {self.amount_of_synapses_in_total}')
            self.logger.info(f'\n  - Participating Nodes are neurons in a tuple with at least:'
                             f' {self.synapse_amount_threshold} synapses')
            self.logger.info(f'  - Synapses in the graph: {self.amount_of_synapses_in_graph}')

        self.logger.info(f'  - Nodes: {len(self.graph)}')
        self.logger.info(f'  - Edges: {len(self.graph.edges)}')
        self.logger.info(f'  - Average clustering coefficient: {round(nx.average_clustering(self.graph), 3)}')

        un_dir_graph = nx.Graph(self.graph)
        avg_short_path_len = np.mean([nx.average_shortest_path_length(un_dir_graph.subgraph(c).copy()) for c in
                                      nx.connected_components(un_dir_graph)])

        self.logger.info(f'  - Average shortest path (undirected): {round(avg_short_path_len, 3)}')
        self.logger.info(f'  - Density: {round(nx.density(self.graph), 3)}')

        max_node, max_degree = sorted(list(self.graph.out_degree), key=lambda x: x[1], reverse=True)[0]
        max_node = self.neuron_names[max_node] if self.neuron_names else max_node
        self.logger.info(f'  - Max out degree of Node {max_node}: {max_degree}')

        max_node, max_degree = sorted(list(self.graph.in_degree), key=lambda x: x[1], reverse=True)[0]
        max_node = self.neuron_names[max_node] if self.neuron_names else max_node
        self.logger.info(f'  - Max in degree of Node {max_node}: {max_degree}')
        self.logger.info('')

        self.plot_properties()

    def node_properties(self, node: Union[str, int]):
        self.logger.info(f'Node {node} properties:')

        if isinstance(node, str):
            node = self.neuron_names.index(node)

        self.logger.info(f'Degree: {self.graph.degree[node]}')
        self.logger.info(f'Out Degree: {self.graph.out_degree[node]}')
        self.logger.info(f'In Degree: {self.graph.in_degree[node]}')
        self.logger.info(f'Clustering coefficient: {round(nx.average_clustering(self.graph, nodes=[node]), 3)}')

    def __plot_rich_club_coefficient(self):
        un_dir_graph = nx.Graph(self.graph)
        un_dir_graph.remove_edges_from(nx.selfloop_edges(un_dir_graph))
        rc = nx.rich_club_coefficient(un_dir_graph, normalized=False, seed=42)

        plt.figure()
        plt.title('Rich Club Coefficient')
        plt.xlabel('Degree (k)')
        plt.ylabel('Rich Club Coefficient')
        plt.scatter(list(rc.keys()), list(rc.values()))
        plt.show()

    def plot_properties(self):
        if not self.plot:
            return
        # TODO: need better plotting tools / motif plotting
        if self.neuron_names:
            mapping = {i: n for i, n in enumerate(self.neuron_names)}
            plot_g = nx.relabel_nodes(self.graph, mapping)
        else:
            plot_g = self.graph

        plt.figure()
        nx.draw_networkx(plot_g, with_labels=True, node_size=600, node_color='lightgreen')
        plt.show()

        self.__plot_rich_club_coefficient()
