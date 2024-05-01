import networkx as nx

from networks.loaders.network_loader import NetworkLoader
from random_networks.markov_chain_switching import MarkovChainSwitching
from utils.config import Config


def graph_generators_comparison():
    loader = NetworkLoader()

    n = 1000
    e = 4000
    p = e / (n * (n - 1))
    er_graph = nx.erdos_renyi_graph(n=n, p=p, directed=True)  # random
    ws_graph = nx.DiGraph(nx.watts_strogatz_graph(n=n, k=10, p=0.5))  # small world, but not scale free
    ba_graph = nx.DiGraph(nx.barabasi_albert_graph(n=n, m=2))  # scale free (undirected), i.e.: power law
    scale_free_graph = nx.DiGraph(nx.scale_free_graph(n=n))  # scale free (directed), i.e.: power law

    print('erdos_renyi')
    loader.load_graph(er_graph)
    print('watts_strogatz')
    loader.load_graph(ws_graph)
    print('barabasi_albert')
    loader.load_graph(ba_graph)
    print('scale_free_graph')
    loader.load_graph(scale_free_graph)


def compare_to_orig_network():
    loader = NetworkLoader()

    network = loader.load_network_file(polarity_xlsx_file_path="networks/data/polarity_2020/s1_data.xlsx",
                                       polarity_sheet_name='5. Sign prediction',
                                       name="polarity 2020 SI 1")

    markov_chain = MarkovChainSwitching(network)
    loader.load_graph(markov_chain.generate(amount=1)[0])

    scale_free_graph = nx.DiGraph(nx.scale_free_graph(n=len(network.graph)), alpha=0.5, beta=0, gamma=0.5, delta_in=0)
    loader.load_graph(scale_free_graph)


if __name__ == "__main__":
    config = Config()
    config.set_property("run_args", "plot_properties", "true")
    config.set_property("run_args", "plot_full_graph", "false")
    config.set_property("run_args", "run_sub_graph_search", "false")
    config.set_property("run_args", "run_motif_criteria", "false")

    # graph_generators_comparison()
    compare_to_orig_network()
