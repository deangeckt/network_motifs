from typing import Union

import networkx as nx
import matplotlib.pyplot as plt
import numpy as np

from post_motif_analysis.polarity_counter import count_network_polarity_ratio
from utils.common import basic_plot
from utils.simple_logger import Logger


class Network:
    def __init__(self, synapse_threshold: int):
        self.logger = Logger()
        self.graph = nx.DiGraph()

        # general configuration
        # TODO: mv plots to notebook
        self.plot_props = False
        self.plot_full_graph = False

        # neurons configuration
        self.neuron_names: list[str] = []
        self.amount_of_synapses_in_graph = 0
        self.amount_of_synapses_in_total = 0
        self.synapse_threshold = synapse_threshold
        self.participating_neurons = set()

        # polarity configuration
        self.use_polarity = False
        self.polarity_ratio = 1  # E/I (+/-) ratio

    def calc_polarity_ratio(self):
        """
        should be called after the graph has been loaded
        """
        polarities = []
        for s, t in self.graph.edges:
            if 'polarity' not in self.graph[s][t]:
                return
            polarities.append(self.graph.get_edge_data(s, t)['polarity'])

        self.polarity_ratio = count_network_polarity_ratio(polarities)
        self.use_polarity = True

    def __degree_stats(self, degree_data: dict, title: str):
        max_node, max_degree = sorted(list(degree_data), key=lambda x: x[1], reverse=True)[0]
        max_node = self.neuron_names[max_node] if self.neuron_names else max_node
        values = np.array(list(dict(degree_data).values()))
        self.logger.info(f'\t{title}: Mean: {round(np.mean(values), 3)} '
                         f'Std: {round(np.std(values), 3)} '
                         f'Median: {round(np.median(values), 3)} '
                         f'Max: {max_degree} (node: {max_node})')

    def properties(self):
        self.logger.info(f'Network properties:')
        if self.neuron_names:
            self.logger.info(f'\tNeurons: {len(self.neuron_names)}')
            self.logger.info(f'\tNeurons with a Synapse: {len(self.participating_neurons)}')
            self.logger.info(f'\tSynapses in the network: {self.amount_of_synapses_in_total}')
            self.logger.info(f'\n\tParticipating Nodes are neurons in a tuple with at least:'
                             f' {self.synapse_threshold} synapses')
            self.logger.info(f'\tSynapses in the graph: {self.amount_of_synapses_in_graph}')

        self.logger.info(f'\tNodes: {len(self.graph)}')
        self.logger.info(f'\tEdges: {len(self.graph.edges)}')
        self.logger.info(f'\tAverage clustering coefficient: {round(nx.average_clustering(self.graph), 3)}')

        un_dir_graph = nx.Graph(self.graph)
        avg_short_path_len = np.mean([nx.average_shortest_path_length(un_dir_graph.subgraph(c).copy()) for c in
                                      nx.connected_components(un_dir_graph)])
        self.logger.info(f'\tAverage shortest path (undirected): {round(avg_short_path_len, 3)}')

        self.logger.info(f'\tDensity: {round(nx.density(self.graph), 3)}')
        if self.use_polarity:
            self.logger.info(f'\tPolarity E/I ratio: {round(self.polarity_ratio, 3)}')

        self.__degree_stats(self.graph.degree, 'Degree')
        self.__degree_stats(self.graph.in_degree, 'In-Degree')
        self.__degree_stats(self.graph.out_degree, 'Out-Degree')

        self.logger.info('')
        self.plot_properties()
        self.plot_graph()

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

    def __plot_degree_dist_log(self, normalized=True):
        y = nx.degree_histogram(self.graph)
        x = np.arange(0, len(y)).tolist()
        n = self.graph.number_of_nodes()

        if normalized:
            for i in range(len(y)):
                y[i] = y[i] / n

        basic_plot(data=(x, y, 'o'),
                   title='Degree Distribution (log-log scale)',
                   xlabel='Degree (log scale)',
                   ylabel='# of Nodes (log scale)',
                   plot_func=plt.plot,
                   log_scale=True)

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
        if not self.plot_props:
            return

        self.__plot_degree_dist()
        self.__plot_degree_dist_log()
        self.__plot_degree_in_dist()
        self.__plot_degree_out_dist()
        # self.__plot_rich_club_coefficient()

    def plot_graph(self):
        if not self.plot_full_graph:
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
