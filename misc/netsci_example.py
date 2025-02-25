import time

import networkx as nx
import numpy as np
import netsci.visualization as nsv
from matplotlib import pyplot as plt
import netsci.metrics.motifs as nsm
from networkx import DiGraph

from utils.export_import import import_network
from utils.types import NetworkBinaryFile


# https://github.com/gialdetti/netsci

# TODO: nd array loader?

def netsci_motif(graph: DiGraph):
    def sanity(A):
        for n_idx, n in enumerate(sorted_nodes):
            if not list(np.array(list(participating_nodes))[np.where(A[n_idx])]) == list(graph.adj[n]):
                print(n)
        print('all good?')

    # correction for cook matrix
    def access_name():
        for n_idx, n in enumerate(sorted_nodes):
            print(f'np idx: {n_idx}, node: {n}, name: {neuron_names[n]}')
    # graph.remove_edges_from(nx.selfloop_edges(graph))

    sorted_nodes = list(graph.nodes)
    sorted_nodes.sort()

    A = nx.adjacency_matrix(graph, nodelist=sorted_nodes).todense()

    # sanity(A)
    # access_name()

    start_time = time.time()
    f, participating = nsm.motifs(A, algorithm='louzoun', participation=True)
    end_time = time.time()
    print(f'timer [Sec]: {round(end_time - start_time, 2)}')

    f_ = f[3:]
    participating = participating[3:]

    motif_keys = [12, 36, 6, 38, 14, 74, 98, 78, 102, 46, 108, 110, 238]
    motifs_dict = {}
    for i, k in enumerate(motif_keys):
        motifs_dict[k] = {'n_real': f_[i], 'sub_graphs': participating[i]}

    motifs_dict = dict(sorted(motifs_dict.items()))
    for k in motifs_dict:
        print(k, motifs_dict[k]['n_real'])
        print(motifs_dict[k]['sub_graphs'])
        print('--')
        for sub_graph in motifs_dict[k]['sub_graphs']:
            sub_graph_neurons = [sorted_nodes[n] for n in sub_graph]
            print(sub_graph_neurons)

        print()

    # nsv.bar_motifs(f)
    # plt.show()


network_binary: NetworkBinaryFile = import_network('../test_network.bin')
neuron_names = network_binary['neuron_names']
participating_nodes = network_binary['participating_nodes']
graph = network_binary['graph']

netsci_motif(graph)
