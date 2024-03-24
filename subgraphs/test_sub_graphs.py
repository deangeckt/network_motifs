from networkx import DiGraph
from tqdm import tqdm

from network import Network
from network_randomizer import NetworkRandomizer
from subgraphs.mfinder_enumeration import MFinder
from subgraphs.specific_subgraphs import SpecificSubGraphs
from subgraphs.sub_graphs_utils import get_sub_id_name, MotifName
from utils.simple_logger import Logger


def __ffl_compare(graph: DiGraph, k: int):
    mfinder = MFinder(graph)
    mfinder_sub_graphs = mfinder.search_sub_graphs(k=k)

    mfinder_ffl_count = 0
    for sub_id in mfinder_sub_graphs:
        sub_name = get_sub_id_name(sub_id=sub_id, k=k)
        if sub_name != 'feed_forwards':
            continue
        mfinder_ffl_count = mfinder_sub_graphs[sub_id]

    specific = SpecificSubGraphs(graph, search=[MotifName.feed_forwards])
    assert mfinder_ffl_count == specific.search_sub_graphs(k=k)[MotifName.feed_forwards]


def test_ffl():
    logger = Logger()
    logger.toggle(False)

    k = 3
    with open("networks/paper_exmp_network.txt", "r") as f:
        network = Network()
        network.load_from_txt(f.readlines())

        print('compare on original network')
        __ffl_compare(network.graph, k)

        rand_amount = 10
        print(f'compare on {rand_amount} random networks')
        randomizer = NetworkRandomizer(network.graph)
        random_networks = randomizer.generate(amount=rand_amount)
        [__ffl_compare(rand_network, k) for rand_network in tqdm(random_networks)]
