import networkx as nx
import pytest

from networks.loaders.network_loader import NetworkLoader
from subgraphs.fanmod_esu import FanmodESU
from subgraphs.mfinder_enum_induced import MFinderInduced
from subgraphs.mfinder_enum_none_induced import MFinderNoneInduced
from isomorphic.isomorphic import get_fsl_ids_iso_mapping, IsomorphicMotifMatch
from utils.sub_graphs import get_sub_id_name, MotifName
from subgraphs.triadic_census import TriadicCensus
from utils.types import SubGraphSearchResult, NetworkInputType, NetworkLoaderArgs

simple_input_args = NetworkLoaderArgs(
    synapse_threshold=5
)


def __compare(k: int, expected_sub_graphs: dict, actual_sub_graphs: SubGraphSearchResult, test_specific=True):
    for sub_id in actual_sub_graphs.fsl:
        if test_specific:

            sub_name = get_sub_id_name(sub_id=sub_id, k=k)
            if sub_name not in expected_sub_graphs:
                continue

            assert actual_sub_graphs.fsl[sub_id] == expected_sub_graphs[sub_name]
        else:
            assert actual_sub_graphs.fsl[sub_id] == expected_sub_graphs[sub_id]


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
    loader = NetworkLoader(simple_input_args)
    network = loader.load_network_file(file_path=network_file,
                                       input_type=NetworkInputType.simple_adj_txt)

    iso_matcher = IsomorphicMotifMatch(k=k, polarity_options=[])
    isomorphic_mapping = iso_matcher.isomorphic_mapping

    mfinder = MFinderInduced(network.graph, isomorphic_mapping)
    mfinder_sub_graphs = mfinder.search_sub_graphs(k=k, allow_self_loops=False)
    __compare(k, expected, mfinder_sub_graphs)

    fanmod = FanmodESU(network.graph, isomorphic_mapping)
    fanmod_sub_graphs = fanmod.search_sub_graphs(k=k, allow_self_loops=False)
    __compare(k, expected, fanmod_sub_graphs)

    triadic_census = TriadicCensus(network.graph, isomorphic_mapping)
    triadic_census_sub_graphs = triadic_census.search_sub_graphs(k=k, allow_self_loops=False)
    __compare(k, expected, triadic_census_sub_graphs)


@pytest.mark.parametrize("network_file,expected", [paper_example_none_induced, paper_ecoli_none_induced])
def test_three_sub_graphs_none_induced(network_file, expected):
    k = 3
    loader = NetworkLoader(simple_input_args)
    network = loader.load_network_file(file_path=network_file,
                                       input_type=NetworkInputType.simple_adj_txt)

    iso_matcher = IsomorphicMotifMatch(k=k, polarity_options=[])
    isomorphic_mapping = iso_matcher.isomorphic_mapping

    mfinder = MFinderNoneInduced(network.graph, isomorphic_mapping)
    mfinder_sub_graphs = mfinder.search_sub_graphs(k=k, allow_self_loops=False)
    __compare(k, expected, mfinder_sub_graphs)


def test_k_2_with_self_loops():
    k = 2

    iso_matcher = IsomorphicMotifMatch(k=k, polarity_options=[], allow_self_loops=True)
    isomorphic_mapping = iso_matcher.isomorphic_mapping

    graph = nx.DiGraph([(1, 2), (1, 1), (1, 3), (3, 2), (3, 4), (4, 4)])
    expected = {2: 1, 3: 2, 5: 1}

    mfinder = MFinderInduced(graph, isomorphic_mapping)
    mfinder_sub_graphs = mfinder.search_sub_graphs(k=k, allow_self_loops=True)
    __compare(k, expected, mfinder_sub_graphs, test_specific=False)

    fanmod = FanmodESU(graph, isomorphic_mapping)
    fanmod_sub_graphs = fanmod.search_sub_graphs(k=k, allow_self_loops=True)
    __compare(k, expected, fanmod_sub_graphs, test_specific=False)


def test_k_2_with_self_loops_wo_iso_mapping():
    k = 2
    graph = nx.DiGraph([(1, 2), (1, 1), (1, 3), (3, 2), (3, 4), (4, 4)])
    expected = {2: 1, 3: 2, 5: 1}

    mfinder = MFinderInduced(graph, {})
    mfinder_sub_graphs = mfinder.search_sub_graphs(k=k, allow_self_loops=True)

    ids_iso_mapping = get_fsl_ids_iso_mapping(list(mfinder_sub_graphs.fsl.keys()), list(expected.keys()), k=k)
    for src_id in ids_iso_mapping:
        tar_id = ids_iso_mapping[src_id]
        assert mfinder_sub_graphs.fsl[src_id] == expected[tar_id]

    fanmod = FanmodESU(graph, {})
    fanmod_sub_graphs = fanmod.search_sub_graphs(k=k, allow_self_loops=True)
    ids_iso_mapping = get_fsl_ids_iso_mapping(list(fanmod_sub_graphs.fsl.keys()), list(expected.keys()), k=k)
    for src_id in ids_iso_mapping:
        tar_id = ids_iso_mapping[src_id]
        assert fanmod_sub_graphs.fsl[src_id] == expected[tar_id]
