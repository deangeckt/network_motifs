from networkx import DiGraph
from tqdm import tqdm

from network import Network
from network_randomizer import NetworkRandomizer
from subgraphs.mfinder_enumeration import MFinder
from subgraphs.specific_subgraphs import SpecificSubGraphs
from subgraphs.sub_graphs_utils import get_sub_id_name, MotifName
from utils.simple_logger import Logger


def __compare(graph: DiGraph, k: int, compare: list[MotifName]):
    mfinder = MFinder(graph)
    mfinder_sub_graphs = mfinder.search_sub_graphs(k=k)

    specific = SpecificSubGraphs(graph, search=compare)
    expected_sub_graphs = specific.search_sub_graphs(k=k)

    for sub_id in mfinder_sub_graphs:
        sub_name = get_sub_id_name(sub_id=sub_id, k=k)
        if sub_name not in compare:
            continue
        assert mfinder_sub_graphs[sub_id] == expected_sub_graphs[sub_name]


def test_three_sub_graphs():
    logger = Logger()
    logger.toggle(False)
    compare = [MotifName.feed_forwards, MotifName.fan_outs, MotifName.cascades]
    k = 3
    with open("networks/paper_exmp_network.txt", "r") as f:
        network = Network()
        network.load_from_txt(f.readlines())

        print('compare on original network')
        __compare(network.graph, k, compare)

        rand_amount = 50
        print(f'compare on {rand_amount} random networks')
        randomizer = NetworkRandomizer(network.graph)
        random_networks = randomizer.generate(amount=rand_amount)
        [__compare(rand_network, k, compare) for rand_network in tqdm(random_networks)]
