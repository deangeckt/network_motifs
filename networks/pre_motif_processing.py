from networks.loaders.network_loader import NetworkLoader
from networks.network import Network
from utils.export_import import export_network
from utils.types import NetworkInputType, NetworkLoaderArgs, NetworkBinaryFile
import networkx as nx
import operator
from functools import reduce

synapse_threshold = 1
synapse_threshold_union = 5


def graph_union(n1: Network, n2: Network) -> Network:
    n1_participating_neurons = [n for i, n in enumerate(n1.neuron_names) if i in n1.participating_nodes]
    n2_participating_neurons = [n for i, n in enumerate(n2.neuron_names) if i in n2.participating_nodes]

    n1_mapping = {n1.neuron_names.index(n): n for n in n1_participating_neurons}
    n2_mapping = {n2.neuron_names.index(n): n for n in n2_participating_neurons}

    n1_mapping_rev = {n: n1.neuron_names.index(n) for n in n1_participating_neurons}
    n2_mapping_rev = {n: n2.neuron_names.index(n) for n in n2_participating_neurons}

    n1_edges = [(n1_mapping[s], n1_mapping[t]) for s, t in n1.graph.edges]
    n2_edges = [(n2_mapping[s], n2_mapping[t]) for s, t in n2.graph.edges]

    union_edges = set(n1_edges).union(set(n2_edges))
    union_neurons = list(set(list(reduce(operator.concat, union_edges))))

    n1_att_edges = {(s, t): att for (s, t, att) in list(n1.graph.edges(data=True))}
    n2_att_edges = {(s, t): att for (s, t, att) in list(n2.graph.edges(data=True))}

    graph = nx.DiGraph()
    for s, t in union_edges:
        s1 = n1_mapping_rev.get(s, -1)
        t1 = n1_mapping_rev.get(t, -1)

        s2 = n2_mapping_rev.get(s, -1)
        t2 = n2_mapping_rev.get(t, -1)

        att1 = n1_att_edges.get((s1, t1), {})
        att2 = n2_att_edges.get((s2, t2), {})

        uni_att = dict(att1)
        uni_att.update(att2)

        if att1 and att2:
            uni_att['synapse'] = att1.get('synapse', 0) + att2.get('synapse', 0)
            uni_att['gap'] = att1.get('gap', 0) + att2.get('gap', 0)
            uni_att['polarity'] = att1.get('polarity') if att1.get('polarity') else att2.get('polarity')

        if uni_att['synapse'] >= synapse_threshold_union:
            graph.add_edge(union_neurons.index(s), union_neurons.index(t), **att1)

    network = Network(synapse_threshold=synapse_threshold_union)
    network.graph = graph
    network.participating_nodes = set(list(range(len(union_neurons))))
    network.neuron_names = union_neurons

    return network


def graph_intersection(n1: Network, n2: Network) -> Network:
    n1_participating_neurons = [n for i, n in enumerate(n1.neuron_names) if i in n1.participating_nodes]
    n2_participating_neurons = [n for i, n in enumerate(n2.neuron_names) if i in n2.participating_nodes]

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
    network.participating_nodes = set(list(range(len(intersection_neurons))))
    network.neuron_names = intersection_neurons

    return network


simple_input_args = NetworkLoaderArgs(
    synapse_threshold=synapse_threshold,
    filter_syn_type='chem',
    filter_sex_type='herm'
)
loader = NetworkLoader(simple_input_args)

n1 = loader.load_network_file(file_path="networks/data/polarity_2020/s1_data.xlsx",
                              input_type=NetworkInputType.polarity_xlsx)

n2 = loader.load_network_file(file_path="networks/data/intersections/ma_and_cook_si2_herm_chem.bin",
                              input_type=NetworkInputType.binary_network)

int_network = graph_union(n1, n2)
int_network.properties()
export_network('networks/data/intersections/pol_m5_OR__ma_AND_cook_si2_herm_chem.bin',
               NetworkBinaryFile(graph=int_network.graph,
                                 participating_nodes=int_network.participating_nodes,
                                 neuron_names=int_network.neuron_names)
               )
