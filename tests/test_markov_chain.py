import networkx as nx

from networks.loaders.network_loader import NetworkLoader
from random_networks.markov_chain_switching import MarkovChainSwitching
from utils.types import NetworkLoaderArgs, NetworkInputType

network_file = "networks/data/Cook_2019/SI 2 Synapse adjacency matrices.xlsx"


def test_markov_chain():
    """
    Test that the randomizer generated 10 networks and that the success rate is greater than 80%
    """
    loader = NetworkLoader(NetworkLoaderArgs(synapse_threshold=5))
    network = loader.load_network_file(file_path=network_file,
                                       input_type=NetworkInputType.worm_wiring_xlsx)

    randomizer = MarkovChainSwitching(network, switch_factor=10)
    amount_of_networks = 10
    random_networks = randomizer.generate(amount_of_networks)

    assert len(random_networks) == amount_of_networks
    success_rate = randomizer.success_switch / randomizer.markov_chain_num_iterations * amount_of_networks
    assert success_rate > 0.8
