import pytest

from network import Network
from subgraphs.mfinder_enum_induced import MFinderInduced
from subgraphs.mfinder_enum_none_induced import MFinderNoneInduced
from subgraphs.sub_graphs_utils import get_sub_id_name, MotifName, generate_isomorphic_k_sub_graphs
from utils.simple_logger import Logger

logger = Logger()
logger.toggle(False)


def __compare(k: int, expected_sub_graphs: dict, actual_sub_graphs: dict):
    for sub_id in actual_sub_graphs:
        sub_name = get_sub_id_name(sub_id=sub_id, k=k)
        if sub_name not in expected_sub_graphs:
            continue
        assert actual_sub_graphs[sub_id] == expected_sub_graphs[sub_name]


paper_example_induced = (r"networks/Uri_Alon_2002/example.txt",
                         {MotifName.cascades: 10, MotifName.fan_outs: 3, MotifName.feed_forwards: 5}
                         )
paper_ecoli_induced = (r"networks/Uri_Alon_2002/coliInterNoAutoRegVec.txt",
                       {MotifName.cascades: 162, MotifName.fan_outs: 226, MotifName.feed_forwards: 40}
                       )
paper_example_none_induced = (r"networks/Uri_Alon_2002/example.txt",
                              {MotifName.cascades: 15, MotifName.fan_outs: 8, MotifName.feed_forwards: 5}
                              )
paper_ecoli_none_induced = (r"networks/Uri_Alon_2002/coliInterNoAutoRegVec.txt",
                            {MotifName.cascades: 202, MotifName.fan_outs: 266, MotifName.feed_forwards: 40}
                            )


@pytest.mark.parametrize("network_file,expected", [paper_example_induced, paper_ecoli_induced])
def test_three_sub_graphs_induced(network_file, expected):
    k = 3
    network = Network()
    network.load_adj_file(network_file)
    isomorphic_mapping, _ = generate_isomorphic_k_sub_graphs(k=k)

    mfinder = MFinderInduced(network.graph, isomorphic_mapping)
    mfinder_sub_graphs = mfinder.search_sub_graphs(k=k)
    __compare(k, expected, mfinder_sub_graphs)


@pytest.mark.parametrize("network_file,expected", [paper_example_none_induced, paper_ecoli_none_induced])
def test_three_sub_graphs_none_induced(network_file, expected):
    k = 3
    network = Network()
    network.load_adj_file(network_file)
    isomorphic_mapping, _ = generate_isomorphic_k_sub_graphs(k=k)

    mfinder = MFinderNoneInduced(network.graph, isomorphic_mapping)
    mfinder_sub_graphs = mfinder.search_sub_graphs(k=k)
    __compare(k, expected, mfinder_sub_graphs)
