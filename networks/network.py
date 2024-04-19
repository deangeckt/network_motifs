from typing import Union

import networkx as nx
import matplotlib.pyplot as plt
import numpy as np

from utils.common import basic_plot
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

        # polarity configuration
        self.use_polarity = False
        self.polarity_ratio = 1  # E/I (+/-) ratio

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

        avg_degree = sum(dict(self.graph.degree).values()) / len(self.graph)
        self.logger.info(f'  - Average degree: {round(avg_degree, 3)}')

        max_node, max_degree = sorted(list(self.graph.out_degree), key=lambda x: x[1], reverse=True)[0]
        max_node = self.neuron_names[max_node] if self.neuron_names else max_node
        self.logger.info(f'  - Max out degree of Node {max_node}: {max_degree}')

        max_node, max_degree = sorted(list(self.graph.in_degree), key=lambda x: x[1], reverse=True)[0]
        max_node = self.neuron_names[max_node] if self.neuron_names else max_node
        self.logger.info(f'  - Max in degree of Node {max_node}: {max_degree}')

        if self.use_polarity:
            self.logger.info(f'  - Polarity E/I ratio: {round(self.polarity_ratio, 3)}')

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
        data = list(rc.keys()), list(rc.values())
        basic_plot(data=data,
                   title='Rich Club Coefficient',
                   xlabel='Degree (k)',
                   ylabel='Rich Club Coefficient',
                   plot_func=plt.scatter)

    def __plot_degree_dist(self):
        degree_sequence = sorted((d for n, d in self.graph.degree()), reverse=True)
        data = np.unique(degree_sequence, return_counts=True)
        basic_plot(data=data,
                   title='Degree Distribution',
                   xlabel='Degree',
                   ylabel='# of Nodes',
                   plot_func=plt.bar)

    def __plot_degree_in_dist(self):
        degree_sequence = sorted((d for n, d in self.graph.in_degree), reverse=True)
        data = np.unique(degree_sequence, return_counts=True)
        basic_plot(data=data,
                   title='Degree In Distribution',
                   xlabel='Degree in',
                   ylabel='# of Nodes',
                   plot_func=plt.bar)

    def __plot_degree_out_dist(self):
        degree_sequence = sorted((d for n, d in self.graph.out_degree), reverse=True)
        data = np.unique(degree_sequence, return_counts=True)
        basic_plot(data=data,
                   title='Degree Out Distribution',
                   xlabel='Degree out',
                   ylabel='# of Nodes',
                   plot_func=plt.bar)

    def plot_properties(self):
        if not self.plot:
            return
        # TODO: need better plotting tools / motif plotting
        if self.neuron_names:
            mapping = {i: n for i, n in enumerate(self.neuron_names)}
            plot_g = nx.relabel_nodes(self.graph, mapping)
        else:
            plot_g = self.graph

        # plt.figure()
        # nx.draw_networkx(plot_g, with_labels=True, node_size=600, node_color='lightgreen')
        # plt.show()

        self.__plot_degree_dist()
        self.__plot_degree_in_dist()
        self.__plot_degree_out_dist()
        self.__plot_rich_club_coefficient()
