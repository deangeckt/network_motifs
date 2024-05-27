import pytest

from large_subgraphs.single_input_moudle import SingleInputModule
from networks.loaders.network_loader import NetworkLoader
from subgraphs.mfinder_enum_induced import MFinderInduced
from subgraphs.sub_graphs_utils import generate_isomorphic_k_sub_graphs
from subgraphs.triadic_census import TriadicCensus
from utils.types import NetworkInputType, NetworkLoaderArgs

simple_input_args = NetworkLoaderArgs(
    synapse_threshold=5
)

paper_example = r"networks/data/Uri_Alon_2002/example.txt"
paper_ecoli = r"networks/data/Uri_Alon_2002/coliInterNoAutoRegVec.txt"


@pytest.mark.parametrize("network_file", [paper_example, paper_ecoli])
def test_fan_out(network_file):
    k = 3
    loader = NetworkLoader(simple_input_args)
    network = loader.load_network_file(file_path=network_file,
                                       input_type=NetworkInputType.simple_adj_txt)

    isomorphic_mapping, _ = generate_isomorphic_k_sub_graphs(k=k)

    triadic_census = TriadicCensus(network.graph, isomorphic_mapping)
    triadic_census_sub_graphs = triadic_census.search_sub_graphs(k=k)

    sim = SingleInputModule(network.graph)
    sim_res = sim.search_sub_graphs(min_control_size=k-1, max_control_size=5)

    assert sim_res.fsl['SIM_2'] == triadic_census_sub_graphs.fsl.get(6, 0)


@pytest.mark.parametrize("network_file", [paper_example, paper_ecoli])
def test_sim3(network_file):
    k = 4
    loader = NetworkLoader(simple_input_args)
    network = loader.load_network_file(file_path=network_file,
                                       input_type=NetworkInputType.simple_adj_txt)

    isomorphic_mapping, _ = generate_isomorphic_k_sub_graphs(k=k)

    mfinder = MFinderInduced(network.graph, isomorphic_mapping)
    mfinder_sub_graphs = mfinder.search_sub_graphs(k=k)

    sim = SingleInputModule(network.graph)
    sim_res = sim.search_sub_graphs(min_control_size=k - 1, max_control_size=5)
    assert sim_res.fsl['SIM_3'] == mfinder_sub_graphs.fsl.get(14, 0)
