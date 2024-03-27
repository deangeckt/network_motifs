import networkx as nx
import numpy as np
import netsci.visualization as nsv
from matplotlib import pyplot as plt
import netsci.metrics.motifs as nsm
from networkx import DiGraph

# https://github.com/gialdetti/netsci


def netsci_motif(graph: DiGraph):
    def graph_to_netsci_array(g: DiGraph) -> np.ndarray:
        nodes = list(g.nodes)
        nodes.sort()
        return nx.adjacency_matrix(graph, nodelist=nodes).todense()

    A = graph_to_netsci_array(graph)
    f = nsm.motifs(A, algorithm='brute-force', participation=False)
    print(f)

    f_ = f
    feed_forward = f_[6]
    cascade = f_[3]
    fanout = f_[5]

    print(f'feed_forward: {feed_forward} cascade: {cascade} fanout: {fanout}')

    nsv.bar_motifs(f)
    plt.show()
