import ast
import random
from typing import Any

import networkx as nx
from networkx import DiGraph
from tabulate import tabulate
from tqdm import tqdm

from motif_criteria import MotifCriteria
from networks.loaders.network_loader import NetworkLoader
from post_motif_analysis.node_counter import sort_node_appearances_in_sub_graph, sort_node_roles_in_sub_graph
from post_motif_analysis.polarity_counter import get_polarity_frequencies, get_all_sub_graph_polarities
from random_networks.erdos_renyi_forced_edges import ErdosRenyiForcedEdges
from random_networks.markov_chain_switching import MarkovChainSwitching
from subgraphs.fanmod_esu import FanmodESU
from subgraphs.mfinder_enum_induced import MFinderInduced
from subgraphs.mfinder_enum_none_induced import MFinderNoneInduced
from subgraphs.specific_subgraphs import SpecificSubGraphs
from subgraphs.sub_graphs_abc import SubGraphsABC
from subgraphs.sub_graphs_utils import generate_isomorphic_k_sub_graphs, create_base_motif
from subgraphs.triadic_census import TriadicCensus
from utils.config import Config
from utils.simple_logger import Logger, LogLvl
import time
import argparse

from utils.types import SubGraphAlgoName, RandomGeneratorAlgoName, NetworkInputType

sub_graph_algorithms = {
    SubGraphAlgoName.specific: SpecificSubGraphs,
    SubGraphAlgoName.mfinder_induced: MFinderInduced,
    SubGraphAlgoName.mfinder_none_induced: MFinderNoneInduced,
    SubGraphAlgoName.fanmod_esu: FanmodESU,
    SubGraphAlgoName.triadic_census: TriadicCensus
}

random_generator_algorithms = {
    RandomGeneratorAlgoName.markov_chain_switching: MarkovChainSwitching,
    RandomGeneratorAlgoName.erdos_renyi: ErdosRenyiForcedEdges
}


def my_input_files():
    # path = "networks/data/polarity_2020/s1_data.xlsx"
    # sheet_name = '5. Sign prediction'
    # file_type = 'polarity_xlsx'
    #
    # path = "networks/data/Cook_2019/SI 2 Synapse adjacency matrices.xlsx"
    # sheet_name = 'herm chem synapse adjacency'  # 'herm gap jn synapse adjacency'
    # file_type = 'worm_wiring_xlsx'
    #
    # path = "networks/data/Cook_2019/SI 5 Connectome adjacency matrices, corrected July 2020.xlsx"
    # sheet_name = 'hermaphrodite chemical'
    # file_type = 'worm_wiring_xlsx'
    #
    # path = "networks/data/Uri_Alon_2002/example.txt"
    # sheet_name = None
    # file_type = 'simple_adj_txt'

    # path = "networks/data/Durbin_1986/neurodata.txt"
    # sheet_name = None
    # file_type = 'durbin_txt'
    pass


def parse_graph(input_: str) -> DiGraph:
    try:
        tuple_list = [tuple(element.replace(' ', '')) for element in input_]
        return nx.DiGraph(tuple_list)
    except:
        raise argparse.ArgumentTypeError(f'parse graph failed, input should be in the format: "1 2", "2 3"')


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("-rs", "--random_seed", help="random seed for the entire program", default=42)

    # [Output files]
    parser.add_argument("-lf", "--log_file", help="file path to save log results", default=None)
    parser.add_argument("-bf", "--bin_file", help="file path to save binary results", default=None)

    # [Input file]
    parser.add_argument("-it", "--input_type",
                        help="file type of the input network",
                        default='graph',
                        choices=['simple_adj_txt', 'worm_wiring_xlsx', 'polarity_xlsx', 'durbin_txt', 'graph'],
                        required=False)
    parser.add_argument("-inf", "--input_network_file",
                        help="file path of the input network",
                        default="networks/data/Durbin_1986/neurodata.txt"
                        )
    parser.add_argument("-ing", "--input_network_graph",
                        help='a graph: list of tuples where each is an edge. in the format: "1 2" "2 3"...',
                        default=['1 2', '1 3'],
                        nargs='+'
                        )
    parser.add_argument("-sn", "--sheet_name",
                        default=None,
                        help="sheet name of an xlsx input network file")

    # [Run args]
    parser.add_argument("-k", "--k",
                        help="the size of sub-graph / motif",
                        default=3)
    parser.add_argument("-sa", "--sub_graph_algorithm",
                        help="sub-graph enumeration algorithm",
                        default='mfinder_i',
                        choices=['mfinder_i', 'mfinder_ni', 'fanmod', 'triadic_census', 'specific'])
    parser.add_argument("-la", "--log_all_motif_sub_graphs",
                        help="log all motifs sub graphs found",
                        default=True)
    # TODO: add bool flag to use iso mapping sub-graphs search
    # implications: without it - you won't get anti-motifs!
    parser.add_argument("-asl", "--allow_self_loops",
                        help="allow self loops in the (pre motif search) isomorphic sub-graphs search",
                        default=False)

    # [Neuronal]
    parser.add_argument("-st", "--synapse_threshold",
                        help="filter neurons with >= # synapses (only in neuron networks files)",
                        default=5)

    # [Polarity]
    parser.add_argument("-fp", "--filter_polarity",
                        help="polarity: filter neurons with polarity",
                        choices=['+', '-', 'no pred', 'complex'],
                        default=['+', '-'],
                        nargs='+')
    parser.add_argument("-fpn", "--filter_prim_nt",
                        help="polarity: filter neurons with primary neurotransmitter",
                        choices=['Glu', 'GABA', 'ACh', 0],
                        default=['Glu', 'GABA', 'ACh', 0],
                        nargs='+')

    # [Durbin]
    parser.add_argument("-dfr", "--durbin_filter_recon",
                        help="durbin: filter animal: JSH: L4 male, N2U: hermaphrodite adult",
                        choices=['N2U', 'N2U'],
                        default='N2U')
    parser.add_argument("-dfs", "--durbin_filter_syn_type",
                        help="durbin: filter synapse type",
                        choices=['chem', 'gap', 'all'],
                        default='chem')

    # [Randomizer]
    parser.add_argument("-r", "--randomizer",
                        help="main randomizer algorithm in a full motif search",
                        default='markov_chain',
                        choices=['markov_chain', 'erdos_renyi'])
    parser.add_argument("-na", "--network_amount",
                        help="amount of random networks to generate in a full motif search",
                        default=1000)
    parser.add_argument("-sf", "--switch_factor",
                        help="number of switch factors done by the markov chain randomizer",
                        default=10)

    # [Motif criteria]
    parser.add_argument("-a", "--alpha",
                        help="motif criteria alpha for testing p value significance",
                        default=0.01)
    parser.add_argument("-ft", "--frequency_threshold",
                        help="motif criteria frequency threshold test",
                        default=0.1)
    parser.add_argument("-ut", "--uniqueness_threshold",
                        help="motif criteria uniqueness threshold test",
                        default=3)
    parser.add_argument("-uut", "--use_uniq_criteria",
                        help="whether to use the uniqueness test",
                        default=False)

    args = parser.parse_args()

    # for arg in vars(args):
    #     print(arg, getattr(args, arg))
    # TODO: log config / args function - nicely

    return args


if __name__ == "__main__":
    args = parse_args()

    random.seed(args.random_seed)
    sub_graph_algo_choice = SubGraphAlgoName(args.sub_graph_algorithm)
    random_generator_algo_choice = RandomGeneratorAlgoName(args.randomizer)

    loader = NetworkLoader()

    input_type = NetworkInputType(args.input_type)

    if input_type == NetworkInputType.graph:
        graph = parse_graph(args.input_network_graph)
        network = loader.load_graph(graph)
    else:
        network = loader.load_network_file(input_type=NetworkInputType(args.input_type),
                                           file_path=args.input_network,
                                           sheet_name=args.sheet_name)
