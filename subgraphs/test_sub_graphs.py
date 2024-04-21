import pytest

from networks.loaders.network_loader import NetworkLoader
from subgraphs.mfinder_enum_induced import MFinderInduced
from subgraphs.mfinder_enum_none_induced import MFinderNoneInduced
from subgraphs.sub_graphs_utils import get_sub_id_name, MotifName, generate_isomorphic_k_sub_graphs
from utils.simple_logger import Logger
from utils.types import SubGraphSearchResult

logger = Logger()
logger.toggle(False)


def __compare(k: int, expected_sub_graphs: dict, actual_sub_graphs: SubGraphSearchResult):
    for sub_id in actual_sub_graphs.fsl:
        sub_name = get_sub_id_name(sub_id=sub_id, k=k)
        if sub_name not in expected_sub_graphs:
            continue
        assert actual_sub_graphs.fsl[sub_id] == expected_sub_graphs[sub_name]


paper_example_induced = (r"networks/data/Uri_Alon_2002/example.txt",
                         {MotifName.cascade: 10, MotifName.fan_out: 3, MotifName.feed_forward: 5}
                         )
paper_ecoli_induced = (r"networks/data/Uri_Alon_2002/coliInterNoAutoRegVec.txt",
                       {MotifName.cascade: 162, MotifName.fan_out: 226, MotifName.feed_forward: 40}
                       )
paper_example_none_induced = (r"networks/data/Uri_Alon_2002/example.txt",
                              {MotifName.cascade: 15, MotifName.fan_out: 8, MotifName.feed_forward: 5}
                              )
paper_ecoli_none_induced = (r"networks/data/Uri_Alon_2002/coliInterNoAutoRegVec.txt",
                            {MotifName.cascade: 202, MotifName.fan_out: 266, MotifName.feed_forward: 40}
                            )


@pytest.mark.parametrize("network_file,expected", [paper_example_induced, paper_ecoli_induced])
def test_three_sub_graphs_induced(network_file, expected):
    k = 3
    loader = NetworkLoader()

    network = loader.load_network_file(name='', adj_file_path=network_file)
    isomorphic_mapping, _ = generate_isomorphic_k_sub_graphs(k=k)

    mfinder = MFinderInduced(network.graph, isomorphic_mapping)
    mfinder_sub_graphs = mfinder.search_sub_graphs(k=k)
    __compare(k, expected, mfinder_sub_graphs)


@pytest.mark.parametrize("network_file,expected", [paper_example_none_induced, paper_ecoli_none_induced])
def test_three_sub_graphs_none_induced(network_file, expected):
    k = 3
    loader = NetworkLoader()
    network = loader.load_network_file(name='', adj_file_path=network_file)
    isomorphic_mapping, _ = generate_isomorphic_k_sub_graphs(k=k)

    mfinder = MFinderNoneInduced(network.graph, isomorphic_mapping)
    mfinder_sub_graphs = mfinder.search_sub_graphs(k=k)
    __compare(k, expected, mfinder_sub_graphs)
