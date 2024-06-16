import random
from typing import Union

import networkx as nx
from networkx import DiGraph
from tqdm import tqdm

from motif_criteria import MotifCriteria
from networks.loaders.network_loader import NetworkLoader
from networks.network import Network
from utils.node_counter import sort_node_appearances_in_sub_graph, sort_node_roles_in_sub_graph
from utils.polarity_counter import get_polarity_frequencies
from random_networks.barabasi_albert_forced_edges import BarabasiAlbertForcedEdges
from random_networks.erdos_renyi_forced_edges import ErdosRenyiForcedEdges
from random_networks.markov_chain_switching import MarkovChainSwitching
from subgraphs.fanmod_esu import FanmodESU
from subgraphs.mfinder_enum_induced import MFinderInduced
from subgraphs.mfinder_enum_none_induced import MFinderNoneInduced
from large_subgraphs.single_input_moudle import SingleInputModule
from subgraphs.specific_subgraphs import SpecificSubGraphs
from subgraphs.sub_graphs_abc import SubGraphsABC
from utils.isomorphic import generate_isomorphic_k_sub_graphs, get_fsl_ids_iso_mapping
from utils.sub_graphs import create_base_motif, create_sim_motif
from subgraphs.triadic_census import TriadicCensus
from utils.export_import import export_results
from utils.logs import log_motif_results, log_sub_graph_args, log_randomizer_args, log_motifs_table
from utils.simple_logger import Logger
import time
import argparse
from argparse import Namespace

from utils.types import SubGraphAlgoName, RandomGeneratorAlgoName, NetworkInputType, NetworkLoaderArgs, \
    MotifCriteriaArgs, Motif, SubGraphSearchResult, BinaryFile, MotifType

sub_graph_algorithms = {
    SubGraphAlgoName.specific: SpecificSubGraphs,
    SubGraphAlgoName.mfinder_induced: MFinderInduced,
    SubGraphAlgoName.mfinder_none_induced: MFinderNoneInduced,
    SubGraphAlgoName.fanmod_esu: FanmodESU,
    SubGraphAlgoName.triadic_census: TriadicCensus
}

random_generator_algorithms = {
    RandomGeneratorAlgoName.markov_chain_switching: MarkovChainSwitching,
    RandomGeneratorAlgoName.erdos_renyi: ErdosRenyiForcedEdges,
    RandomGeneratorAlgoName.barabasi: BarabasiAlbertForcedEdges
}


def my_input_files():
    # path = "networks/data/Cook_2019/SI 2 Synapse adjacency matrices.xlsx"
    # path = "networks/data/Cook_2019/SI 5 Connectome adjacency matrices, corrected July 2020.xlsx"
    # path = "networks/data/polarity_2020/s1_data.xlsx"
    # path = "networks/data/Uri_Alon_2002/example.txt"
    # path = "networks/data/Durbin_1986/neurodata.txt"
    pass


def parse_graph(input_: str) -> DiGraph:
    tuple_list = [tuple(element.replace(' ', '')) for element in input_]
    return nx.DiGraph(tuple_list)


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("-rs", "--random_seed", help="random seed for the entire program", default=42)
    # [Output files]
    parser.add_argument("-lf", "--log_file",
                        help="file path to save log results",
                        default=None)
    parser.add_argument("-bf", "--bin_file",
                        help="file path to save binary results",
                        default="results/pol_k3_m25.bin")

    # [Input file]
    parser.add_argument("-it", "--input_type",
                        help="the type of the input network",
                        default='polarity_xlsx',
                        choices=['simple_adj_txt', 'worm_wiring_xlsx', 'polarity_xlsx', 'durbin_txt', 'graph'],
                        required=False)
    parser.add_argument("-inf", "--input_network_file",
                        help="file path of the input network",
                        default="networks/data/polarity_2020/s1_data.xlsx"
                        )
    parser.add_argument("-ing", "--input_network_graph",
                        help='a graph: list of tuples where each is an edge. in the format: "1 2" "2 3"...',
                        default=None,
                        nargs='+'
                        )

    # [Run args]
    parser.add_argument("-rmc", "--run_motif_criteria",
                        help="run full motif search with motif criteria tests",
                        action='store_true',
                        default=True)
    parser.add_argument("-sa", "--sub_graph_algorithm",
                        help="sub-graph enumeration algorithm",
                        default='mfinder_i',
                        choices=['mfinder_i', 'mfinder_ni', 'fanmod', 'triadic_census', 'specific'])
    parser.add_argument("-k", "--k",
                        help="the size of sub-graph / motif to search in the enumeration algorithm",
                        type=int,
                        default=3)
    parser.add_argument("-sim", "--sim",
                        help="the maximum size of control size in the SIM search algorithm",
                        type=int,
                        default=3)

    parser.add_argument("-uim", "--use_isomorphic_mapping",
                        help="run (pre motif search) isomorphic sub-graphs search",
                        action='store_true',
                        default=True)
    parser.add_argument("-asl", "--allow_self_loops",
                        help="allow self loops in the (pre motif search) isomorphic sub-graphs search",
                        action='store_true',
                        default=False)

    # [Neuronal]
    parser.add_argument("-st", "--synapse_threshold",
                        help="filter neurons with >= # synapses (only in neuron networks files)",
                        type=int,
                        default=25)
    parser.add_argument("-fsy", "--filter_syn_type",
                        help="filter synapse type, supported in durbin and worm_wiring networks",
                        choices=['chem', 'gap', 'all'],
                        default='chem')
    parser.add_argument("-fsx", "--filter_sex_type",
                        help="filter sex type, supported in durbin and worm_wiring networks",
                        choices=['herm', 'male'],
                        default='herm')

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

    # [Randomizer]
    parser.add_argument("-r", "--randomizer",
                        help="main randomizer algorithm in a full motif search",
                        default='markov_chain',
                        choices=['markov_chain', 'erdos_renyi', 'barabasi'])
    parser.add_argument("-na", "--network_amount",
                        help="amount of random networks to generate in a full motif search",
                        type=int,
                        default=1000)
    parser.add_argument("-sf", "--switch_factor",
                        help="number of switch factors done by the markov chain randomizer",
                        type=int,
                        default=10)

    # [Motif criteria]
    parser.add_argument("-a", "--alpha",
                        help="motif criteria alpha for testing p value significance",
                        type=float,
                        default=0.01)
    parser.add_argument("-ft", "--frequency_threshold",
                        help="motif criteria frequency threshold test",
                        type=float,
                        default=0.1)
    parser.add_argument("-ut", "--uniqueness_threshold",
                        help="motif criteria uniqueness threshold test",
                        type=int,
                        default=3)
    parser.add_argument("-uut", "--use_uniq_criteria",
                        help="whether to use the uniqueness test",
                        action='store_true',
                        default=False)

    return parser.parse_args()


def polarity_motif_search(
        motif_candidates: dict[int, Motif],
        random_network_sub_graph_results: list[SubGraphSearchResult]):
    if not network.use_polarity:
        return

    print('starting pol motif search:')
    for sub_id in tqdm(motif_candidates):
        motif = motif_candidates[sub_id]

        # count polarity frequencies for the random networks
        random_network_polarity_frequencies = []
        for rand_network_res in random_network_sub_graph_results:
            random_network_polarity_frequencies.append(
                get_polarity_frequencies(appearances=rand_network_res.fsl_fully_mapped.get(sub_id, []),
                                         roles=motif.role_pattern,
                                         polarity_options=network.polarity_options,
                                         motif_id=sub_id
                                         ))

        for polarity_motif in motif.polarity_motifs:
            random_network_samples: list[int] = []
            for rand_network_pol_freq in random_network_polarity_frequencies:
                for rand_pol_freq in rand_network_pol_freq:
                    if rand_pol_freq.polarity == polarity_motif.polarity:
                        random_network_samples.append(rand_pol_freq.frequency)
                        break

            polarity_motif.random_network_samples = random_network_samples
            polarity_motif.motif_criteria = motif_criteria.is_motif(polarity_motif)

        log_motifs_table([m for m in motif.polarity_motifs if m.motif_criteria.is_motif != MotifType.none])

    if args.bin_file:
        export_results(BinaryFile(args=args, motifs=motif_candidates))


def _populate_motif(motif: Motif, sub_graphs: list):
    motif.node_roles = sort_node_roles_in_sub_graph(appearances=sub_graphs,
                                                    neuron_names=network.neuron_names,
                                                    motif=motif,
                                                    )
    motif.node_appearances = sort_node_appearances_in_sub_graph(appearances=sub_graphs,
                                                                neuron_names=network.neuron_names)


def sub_graph_search(args: Namespace) -> dict[Union[str, int], Motif]:
    log_sub_graph_args(args)

    sub_graph_algo: SubGraphsABC = sub_graph_algorithms[sub_graph_algo_choice](network.graph, isomorphic_mapping)
    start_time = time.time()
    search_result = sub_graph_algo.search_sub_graphs(k=args.k, allow_self_loops=args.allow_self_loops)
    end_time = time.time()
    logger.info(f'Sub Graph search timer [Sec]: {round(end_time - start_time, 2)}')

    sim = SingleInputModule(network.graph)
    start_time = time.time()
    sim_search_result = sim.search_sub_graphs(min_control_size=args.k, max_control_size=args.sim)
    end_time = time.time()
    logger.info(f'SIM search timer [Sec]: {round(end_time - start_time, 2)}')

    motifs = {}
    loop_over = isomorphic_graphs if isomorphic_graphs else search_result.fsl
    for sub_id in loop_over:
        motif = create_base_motif(sub_id=sub_id, k=args.k)
        motif.n_real = search_result.fsl.get(sub_id, 0)
        motif.sub_graphs = search_result.fsl_fully_mapped.get(sub_id, [])
        _populate_motif(motif=motif, sub_graphs=motif.sub_graphs)
        motifs[sub_id] = motif

    for sim_id in sim_search_result.fsl:
        motif = create_sim_motif(sim_id=sim_id, adj_mat=sim_search_result.adj_mat[sim_id])
        motif.n_real = sim_search_result.fsl[sim_id]
        motif.sub_graphs = sim_search_result.fsl_fully_mapped[sim_id]
        _populate_motif(motif=motif, sub_graphs=motif.sub_graphs)
        motifs[sim_id] = motif

    start_time = time.time()
    if network.use_polarity:
        for sub_id in motifs:
            motif = motifs[sub_id]
            polarity_frequencies = get_polarity_frequencies(appearances=motif.sub_graphs,
                                                            roles=motif.role_pattern,
                                                            polarity_options=network.polarity_options,
                                                            motif_id=sub_id)
            for motif_pol_freq in polarity_frequencies:
                # TODO: compare to 'sim', 'dor' etc... remove isinstance. and add SIM's super class...
                if isinstance(sub_id, str):
                    polarity_motif = create_sim_motif(sim_id=sub_id, adj_mat=sim_search_result.adj_mat[sub_id])
                else:
                    polarity_motif = create_base_motif(sub_id=sub_id, k=args.k)

                polarity_motif.polarity = motif_pol_freq.polarity
                polarity_motif.id = f'{sub_id} {str(motif_pol_freq.polarity)}'
                polarity_motif.n_real = motif_pol_freq.frequency
                polarity_motif.sub_graphs = motif_pol_freq.sub_graphs
                _populate_motif(motif=polarity_motif, sub_graphs=motif_pol_freq.sub_graphs)
                motif.polarity_motifs.append(polarity_motif)
    end_time = time.time()
    logger.info(f'Populate polarity motifs timer [Sec]: {round(end_time - start_time, 2)}')

    return motifs


def motif_search(args: Namespace):
    motif_candidates = sub_graph_search(args)

    if not args.run_motif_criteria:
        log_motif_results(motif_candidates)
        if args.bin_file:
            export_results(BinaryFile(args=args, motifs=motif_candidates))
        return

    log_randomizer_args(args)
    random_generator_algo_choice = RandomGeneratorAlgoName(args.randomizer)
    if random_generator_algo_choice == RandomGeneratorAlgoName.markov_chain_switching:
        randomizer = MarkovChainSwitching(network, switch_factor=args.switch_factor)
    else:
        randomizer = random_generator_algorithms[random_generator_algo_choice](network)

    random_network_amount = args.network_amount
    random_networks = randomizer.generate(amount=random_network_amount)

    random_network_sub_graph_results = []
    for rand_network in tqdm(random_networks):
        sub_graph_algo: SubGraphsABC = sub_graph_algorithms[sub_graph_algo_choice](rand_network, isomorphic_mapping)
        sub_graph_search_result = sub_graph_algo.search_sub_graphs(k=args.k, allow_self_loops=args.allow_self_loops)

        sim = SingleInputModule(rand_network)
        sim_search_result = sim.search_sub_graphs(min_control_size=args.k, max_control_size=args.sim)

        combined_res = SubGraphSearchResult(fsl={**sub_graph_search_result.fsl, **sim_search_result.fsl},
                                            fsl_fully_mapped={**sub_graph_search_result.fsl_fully_mapped,
                                                              **sim_search_result.fsl_fully_mapped})

        random_network_sub_graph_results.append(combined_res)

    # handle mapping in case no isomorphic_graphs
    if not isomorphic_graphs:
        n_real_ids = [m.id for m in motif_candidates.values() if isinstance(m.id, int)]
        rand_network_ids = []
        for rand_network in random_network_sub_graph_results:
            rand_network_ids.append([k for k in list(rand_network.fsl.keys()) if isinstance(k, int)])

        iso_mappings = [get_fsl_ids_iso_mapping(n_real_ids, rand_ids, k=args.k) for rand_ids in rand_network_ids]
        for iso_map, rand_network in zip(iso_mappings, random_network_sub_graph_results):
            for src_ in iso_map:
                tar_ = iso_map[src_]
                if src_ == tar_:
                    continue
                rand_network.fsl[src_] = rand_network.fsl.pop(tar_, 0)
                rand_network.fsl_fully_mapped[src_] = rand_network.fsl_fully_mapped.pop(tar_, [])

    for sub_id in motif_candidates:
        random_network_samples = [rand_network.fsl.get(sub_id, 0) for rand_network in random_network_sub_graph_results]
        motif_candidate: Motif = motif_candidates[sub_id]
        motif_candidate.random_network_samples = random_network_samples
        motif_candidate.motif_criteria = motif_criteria.is_motif(motif_candidate)
        motif_candidates[sub_id] = motif_candidate

    log_motif_results(motif_candidates)

    if args.bin_file and not network.use_polarity:
        export_results(BinaryFile(args=args, motifs=motif_candidates))

    polarity_motif_search(motif_candidates, random_network_sub_graph_results)


def load_network_from_args(args: Namespace) -> Network:
    loader = NetworkLoader(args=NetworkLoaderArgs(**vars(args)))

    input_type = NetworkInputType(args.input_type)
    if input_type == NetworkInputType.graph:
        graph = parse_graph(args.input_network_graph)
        network = loader.load_graph(graph)
    else:
        network = loader.load_network_file(input_type=input_type,
                                           file_path=args.input_network_file)
    return network


if __name__ == "__main__":
    args = parse_args()

    logger = Logger()
    if args.log_file:
        logger.change_file(args.log_file)

    random.seed(args.random_seed)
    sub_graph_algo_choice = SubGraphAlgoName(args.sub_graph_algorithm)

    network = load_network_from_args(args)
    if not len(network.graph):
        logger.info('Empty graph')
        exit(0)

    motif_criteria = MotifCriteria(MotifCriteriaArgs(**vars(args)))
    if args.use_isomorphic_mapping:
        isomorphic_mapping, isomorphic_graphs = generate_isomorphic_k_sub_graphs(k=args.k)
    else:
        isomorphic_mapping, isomorphic_graphs = {}, {}
    motif_search(args)
