from networks.loaders.network_loader import NetworkLoader
from networks.network import Network
from utils.export_import import export_network
from utils.types import NetworkInputType, NetworkLoaderArgs, NetworkBinaryFile
import networkx as nx
import operator
from functools import reduce

synapse_threshold = 1


def graph_intersection(n1: Network, n2: Network) -> Network:
    n1_participating_neurons = [n for i, n in enumerate(n1.neuron_names) if i in n1.participating_neurons]
    n2_participating_neurons = [n for i, n in enumerate(n2.neuron_names) if i in n2.participating_neurons]

    n1_mapping = {n1.neuron_names.index(n): n for n in n1_participating_neurons}
    n2_mapping = {n2.neuron_names.index(n): n for n in n2_participating_neurons}

    n1_mapping_rev = {n: n1.neuron_names.index(n) for n in n1_participating_neurons}
    n2_mapping_rev = {n: n2.neuron_names.index(n) for n in n2_participating_neurons}

    n1_edges = [(n1_mapping[s], n1_mapping[t]) for s, t in n1.graph.edges]
    n2_edges = [(n2_mapping[s], n2_mapping[t]) for s, t in n2.graph.edges]

    intersection_edges = set(n1_edges).intersection(set(n2_edges))
    intersection_neurons = list(set(list(reduce(operator.concat, intersection_edges))))

    n1_att_edges = {(s, t): att for (s, t, att) in list(n1.graph.edges(data=True))}
    n2_att_edges = {(s, t): att for (s, t, att) in list(n2.graph.edges(data=True))}

    graph = nx.DiGraph()
    for s, t in intersection_edges:
        s1 = n1_mapping_rev[s]
        t1 = n1_mapping_rev[t]

        s2 = n2_mapping_rev[s]
        t2 = n2_mapping_rev[t]

        att1 = n1_att_edges[(s1, t1)]
        att2 = n2_att_edges[(s2, t2)]

        att1.update(att2)
        graph.add_edge(intersection_neurons.index(s), intersection_neurons.index(t), **att1)

    network = Network(synapse_threshold=synapse_threshold)
    network.graph = graph
    network.participating_neurons = set(intersection_neurons)
    network.neuron_names = intersection_neurons

    return network


simple_input_args = NetworkLoaderArgs(
    synapse_threshold=synapse_threshold,
    filter_syn_type='chem',
    filter_sex_type='herm'
)
loader = NetworkLoader(simple_input_args)

n1 = loader.load_network_file(file_path="networks/data/Cook_2019/SI 2 Synapse adjacency matrices.xlsx",
                              input_type=NetworkInputType.worm_wiring_xlsx)

n2 = loader.load_network_file(file_path="networks/data/Multilayer_Connectome_2016/edgelist_MA.csv",
                              input_type=NetworkInputType.multilayer)

int_network = graph_intersection(n1, n2)
int_network.properties()
export_network('networks/data/intersections/ma_and_cook_si2_herm_chem.bin',
               NetworkBinaryFile(graph=int_network.graph,
                                 participating_neurons=int_network.participating_neurons,
                                 neuron_names=int_network.neuron_names)
               )
