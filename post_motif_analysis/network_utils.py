from typing import Union

import networkx as nx
import matplotlib.pyplot as plt

from networks.network import Network
import numpy as np


def draw_sub_graph(network: Network, neurons: list[str], full_pol=True, center=None):
    node_list = [network.neuron_names.index(neuron) for neuron in neurons]

    induced_sub_graph = nx.induced_subgraph(network.graph, node_list)
    mapping = {i: n for i, n in enumerate(network.neuron_names)}
    graph = nx.relabel_nodes(induced_sub_graph, mapping)

    if len(graph.edges) == 0:
        print('empty graph')
        return

    pos = nx.circular_layout(graph)
    if center:
        pos[center] = np.array([0, 0])

    full_polarity = nx.get_edge_attributes(graph, 'full_polarity')
    polarity = nx.get_edge_attributes(graph, 'polarity')
    for label in polarity:
        if polarity[label] == 'complex':
            polarity[label] = 'c'

    if full_pol:
        polarity_data = full_polarity
        font_size = 8
    else:
        polarity_data = polarity
        font_size = 10

    synapse = nx.get_edge_attributes(graph, 'synapse')
    gap = nx.get_edge_attributes(graph, 'gap')

    widths = {}
    colors = []
    for edge in synapse:
        widths[edge] = synapse[edge] + gap[edge]
        if synapse[edge] > 0 and gap[edge] > 0:
            colors.append('black')
        elif synapse[edge] > 0:
            colors.append('blue')
        elif gap[edge] > 0:
            colors.append('red')
        else:
            print('err')

    width = list(widths.values())
    width_max = max(width)
    width = list(np.array(width) / width_max)

    # print(synapse)
    # print()
    # print(gap)

    fig, ax = plt.subplots(nrows=1, ncols=1)
    ax.axis('off')

    nx.draw_networkx(graph, pos=pos, ax=ax, node_size=600, node_color='lightgreen', width=width, edge_color=colors,
                     connectionstyle="arc3,rad=0.1")

    nx.draw_networkx_edge_labels(graph, pos=pos, ax=ax, edge_labels=polarity_data, font_color='k', label_pos=0.3,
                                 font_size=font_size)


def draw_neighbors(network: Network, neuron: str, n_type: str, full_pol=True):
    in_neighbors = [network.neuron_names[i] for i in
                    list(network.graph.predecessors(network.neuron_names.index(neuron)))]
    out_neighbors = [network.neuron_names[i] for i in
                     list(network.graph.successors(network.neuron_names.index(neuron)))]

    if n_type == 'in':
        neighbors = in_neighbors
    elif n_type == 'out':
        neighbors = out_neighbors
    else:
        neighbors = in_neighbors
        neighbors.extend(out_neighbors)

    print(neighbors)
    neighbors.append(neuron)
    draw_sub_graph(network, neighbors, full_pol=full_pol, center=neuron)


def node_properties(network: Network, node: Union[str, int]):
    print(f'Node {node} properties:')

    node_idx = network.neuron_names.index(node) if isinstance(node, str) else node
    print(f'Degree: {network.graph.degree[node_idx]}')
    print(f'Out Degree: {network.graph.out_degree[node_idx]}')
    print(f'In Degree: {network.graph.in_degree[node_idx]}')
    print(f'Clustering coefficient: {round(nx.average_clustering(network.graph, nodes=[node_idx]), 3)}')
