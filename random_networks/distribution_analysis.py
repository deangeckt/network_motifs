import networkx as nx

from networks.loaders.network_loader import NetworkLoader
from random_networks.barabasi_albert_forced_edges import BarabasiAlbertForcedEdges
from random_networks.markov_chain_switching import MarkovChainSwitching
from utils.types import NetworkInputType, NetworkLoaderArgs

# TODO: mv to notebook
simple_input_args = NetworkLoaderArgs(
    synapse_threshold=5,
    filter_polarity=['+', '-'],
    filter_prim_nt=['GABA', 'Glu', 'ACh', 0]
)


def graph_generators_comparison():
    loader = NetworkLoader(simple_input_args)

    n = 1000
    e = 4000
    p = e / (n * (n - 1))
    er_graph = nx.erdos_renyi_graph(n=n, p=p, directed=True)  # random

    # small world (high clustering coefficient), but not scale free (power law) i.e: no hubs.
    # undirected. with |E| = n * k/2
    ws_graph = nx.DiGraph(nx.watts_strogatz_graph(n=n, k=10, p=0))

    # small world + scale free but undirected. with |E| = n * m.
    ba_graph = nx.DiGraph(nx.barabasi_albert_graph(n=n, m=2))

    # small world + scale free and directed, but can't control |E|
    scale_free_graph = nx.DiGraph(nx.scale_free_graph(n=n))

    print('erdos_renyi')
    loader.load_graph(er_graph)
    print('watts_strogatz')
    loader.load_graph(ws_graph)
    print('barabasi_albert')
    loader.load_graph(ba_graph)
    print('scale_free_graph')
    loader.load_graph(scale_free_graph)


def compare_to_orig_network():
    loader = NetworkLoader(simple_input_args)
    network = loader.load_network_file(file_path="networks/data/polarity_2020/s1_data.xlsx",
                                       sheet_name='5. Sign prediction',
                                       input_type=NetworkInputType.polarity_xlsx)

    markov_chain = MarkovChainSwitching(network, switch_factor=10)
    loader.load_graph(markov_chain.generate(amount=1)[0])

    barabasi_albert = BarabasiAlbertForcedEdges(network)
    loader.load_graph(barabasi_albert.generate(amount=1)[0])


if __name__ == "__main__":
    # graph_generators_comparison()
    compare_to_orig_network()
