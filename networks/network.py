import collections

import networkx as nx
import numpy as np

from utils.polarity_counter import count_network_polarity_ratio
from utils.simple_logger import Logger


class Network:
    def __init__(self, synapse_threshold: int):
        self.logger = Logger()
        self.graph = nx.DiGraph()

        # neurons configuration
        self.neuron_names: list[str] = []
        self.amount_of_synapses_in_graph = 0
        self.amount_of_synapses_in_total = 0
        self.synapse_threshold = synapse_threshold

        self.amount_of_gaps_in_graph = 0
        self.amount_of_gaps_in_total = 0

        self.participating_neurons = set()

        # polarity configuration
        self.use_polarity = False
        self.polarity_ratio = collections.Counter([])
        self.polarity_options: list[str] = []

        # multilayer configuration
        self.use_monoamines = False

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
        self.polarity_options = list(self.polarity_ratio.keys())
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
        self.logger.info(f'\nNetwork properties:')
        if self.neuron_names:
            self.logger.info(f'\tNeurons in the network: {len(self.neuron_names)}')
            self.logger.info(f'\tParticipating Neurons (in the graph): {len(self.participating_neurons)}')

            self.logger.info(f'\n\tParticipating Nodes are neurons with at least: {self.synapse_threshold} synapses')

            self.logger.info(f'\tSynapses in the network: {self.amount_of_synapses_in_total}')
            self.logger.info(f'\tSynapses in the graph: {self.amount_of_synapses_in_graph}')

            self.logger.info(f'\tGaps in the network: {self.amount_of_gaps_in_total}')
            self.logger.info(f'\tGaps in the graph: {self.amount_of_gaps_in_graph}')

        self.logger.info(f'\tNodes: {len(self.graph)}')
        self.logger.info(f'\tEdges: {len(self.graph.edges)}')

        if not len(self.graph):
            return

        self.logger.info(f'\tAverage clustering coefficient: {round(nx.average_clustering(self.graph), 3)}')

        un_dir_graph = nx.Graph(self.graph)
        avg_short_path_len = np.mean([nx.average_shortest_path_length(un_dir_graph.subgraph(c).copy()) for c in
                                      nx.connected_components(un_dir_graph)])
        self.logger.info(f'\tAverage shortest path (undirected): {round(avg_short_path_len, 3)}')

        self.logger.info(f'\tDensity: {round(nx.density(self.graph), 3)}')
        if self.use_polarity:
            self.logger.info(f'\tPolarity ratios: {self.polarity_ratio}')

        self.__degree_stats(self.graph.degree, 'Degree')
        self.__degree_stats(self.graph.in_degree, 'In-Degree')
        self.__degree_stats(self.graph.out_degree, 'Out-Degree')

    def export_to_gephi(self, output_file_path: str):
        mapping = {i: n for i, n in enumerate(self.neuron_names)}
        g = nx.relabel_nodes(self.graph, mapping) if mapping else self.graph
        nx.write_gexf(g, output_file_path)
