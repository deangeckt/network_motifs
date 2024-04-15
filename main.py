import random

import networkx as nx
from networkx import DiGraph
from tabulate import tabulate
from tqdm import tqdm

from motif_criteria import MotifCriteria
from networks.loaders.network_loader import NetworkLoader
from networks.network import Network
from post_motif_analysis.node_counter import sort_node_appearances_in_sub_graph, sort_node_roles_in_sub_graph
from post_motif_analysis.polarity_counter import get_polarity_frequencies, get_all_sub_graph_polarities
from random_networks.markov_chain_switching import MarkovChainSwitching
from subgraphs.mfinder_enum_induced import MFinderInduced
from subgraphs.mfinder_enum_none_induced import MFinderNoneInduced
from subgraphs.specific_subgraphs import SpecificSubGraphs
from subgraphs.sub_graphs_abc import SubGraphsABC
from subgraphs.sub_graphs_utils import generate_isomorphic_k_sub_graphs, create_base_motif
from utils.config import Config
from utils.simple_logger import Logger, LogLvl
import time

from utils.types import SubGraphAlgoName, Motif, MotifCriteriaResults

sub_graph_algorithms = {
    SubGraphAlgoName.specific: SpecificSubGraphs,
    SubGraphAlgoName.mfinder_induced: MFinderInduced,
    SubGraphAlgoName.mfinder_none_induced: MFinderNoneInduced
}


def log_motif_results(motifs: dict[int, Motif]):
    table = []
    for motif_id in motifs:
        motif = motifs[motif_id]
        motif_criteria = motif.motif_criteria if motif.motif_criteria is not None else MotifCriteriaResults(
            n_real=motif.n_real,
        )

        table.append([motif.id,
                      motif.name,
                      str(motif.adj_mat),
                      motif_criteria.is_motif,
                      motif_criteria.n_real,
                      motif_criteria.n_rand,
                      motif_criteria.std,
                      motif_criteria.z_score,
                      motif_criteria.uniq,
                      motif_criteria.is_statistically_significant,
                      motif_criteria.is_motif_frequent,
                      motif_criteria.is_uniq,
                      motif_criteria.is_anti_motif_frequent
                      ])
    headers = ['ID', 'Name', 'Adj_mat', 'is-motif', 'N_real', 'N_rand', 'Std', 'Z_score', 'uniq',
               'is-significant', 'is-freq', 'is-uniq', 'is-anti-freq']

    col_align = tuple(['center'] * len(headers))
    float_fmt = tuple([".2f"] * len(headers))
    logger.info(tabulate(table, tablefmt="grid", headers=headers, colalign=col_align, floatfmt=float_fmt))

    if not config.get_boolean_property('run_args', 'log_all_motif_sub_graphs'):
        return

    for motif_id in motifs:
        motif = motifs[motif_id]
        if motif.n_real == 0:
            continue
        logger.info(f'\nMotif Id: {motif_id}')

        if not network.use_polarity:
            logger.info(f'Appearances (all sub-graphs): {motif.sub_graphs}')
        else:
            logger.info(
                f'Appearances (all sub-graphs): {get_all_sub_graph_polarities(motif.sub_graphs, network.graph)}')

        logger.info(f'Appearances of Nodes in the sub-graph: {motif.node_appearances}')
        logger.info(f'Role pattern: {motif.role_pattern}')
        for role in motif.node_roles:
            logger.info(f'Appearances of Nodes in role {role}: {motif.node_roles[role]}')
        for pol_freq in motif.polarity_frequencies:
            logger.info(f'Polarity: {pol_freq.polarity} - frequency: {pol_freq.frequency}')


def sub_graph_search(network: Network) -> dict[int, Motif]:
    sub_graph_algo: SubGraphsABC = sub_graph_algorithms[algo](network.graph, isomorphic_mapping)

    logger.info(f'Sub Graph search using Algorithm: {algo}')
    logger.info(f'Sub Graph search using k: {k}')
    logger.info(f'Amount of isomorphic sub graphs of size k={k} is: {len(isomorphic_graphs)}')
    allow_self_loops = config.get_boolean_property('run_args', 'allow_self_loops')
    logger.info(f'Allow self loops: {allow_self_loops}')

    start_time = time.time()
    search_result = sub_graph_algo.search_sub_graphs(k=k)
    end_time = time.time()
    logger.info(f'\nSub Graph search timer [Sec]: {round(end_time - start_time, 2)}')

    total_sub_graphs = sum(search_result.fsl.values())
    logger.info(f'Motif candidates found: {len(search_result.fsl)}')
    logger.info(f'Total number of {k}-node sub graphs found: {total_sub_graphs}')

    motifs = {}
    for sub_id in search_result.fsl:
        motif = create_base_motif(sub_id=sub_id, k=k)
        motif.n_real = search_result.fsl[sub_id]
        motif.sub_graphs = search_result.fsl_fully_mapped[sub_id]
        motif.node_roles = sort_node_roles_in_sub_graph(appearances=motif.sub_graphs,
                                                        neuron_names=network.neuron_names,
                                                        roles=motif.role_pattern)
        motif.node_appearances = sort_node_appearances_in_sub_graph(appearances=motif.sub_graphs,
                                                                    neuron_names=network.neuron_names)
        if network.use_polarity:
            motif.polarity_frequencies = get_polarity_frequencies(appearances=motif.sub_graphs,
                                                                  roles=motif.role_pattern,
                                                                  graph=network.graph)
        motifs[sub_id] = motif

    return motifs


def motif_search(network: Network):
    if not config.get_boolean_property('run_args', 'run_sub_graph_search'):
        return

    motif_candidates = sub_graph_search(network)

    if not config.get_boolean_property('run_args', 'run_motif_criteria'):
        log_motif_results(motif_candidates)
        return

    randomizer = MarkovChainSwitching(network.graph)
    random_network_amount = int(config.get_property('random', 'network_amount'))
    random_networks = randomizer.generate(amount=random_network_amount)

    random_network_sub_graphs = []
    for rand_network in tqdm(random_networks):
        sub_graph_algo: SubGraphsABC = sub_graph_algorithms[algo](rand_network, isomorphic_mapping)
        random_network_sub_graphs.append(sub_graph_algo.search_sub_graphs(k=k))

    motif_criteria = MotifCriteria()
    for sub_id in tqdm(isomorphic_graphs):
        random_network_samples = [rand_network.fsl.get(sub_id, 0) for rand_network in random_network_sub_graphs]
        motif_candidate: Motif = motif_candidates.get(sub_id, create_base_motif(sub_id=sub_id, k=k))
        motif_candidate.random_network_samples = random_network_samples
        motif_candidate.motif_criteria = motif_criteria.is_motif(motif_candidate)
        motif_candidates[sub_id] = motif_candidate

    log_motif_results(motif_candidates)


if __name__ == "__main__":
    random.seed(42)
    logger = Logger(LogLvl.info)
    config = Config()
    loader = NetworkLoader()

    k = int(config.get_property('run_args', 'k'))
    algo = SubGraphAlgoName(config.get_property('run_args', 'sub_graph_algorithm'))
    isomorphic_mapping, isomorphic_graphs = generate_isomorphic_k_sub_graphs(k=k)

    # network = loader.load_graph(nx.DiGraph([(1, 0), (2, 0), (1, 2)]))

    # network = loader.load_network_file(
    #     worm_wiring_xlsx_file_path="networks/data/Cook_2019/SI 2 Synapse adjacency matrices.xlsx",
    #     worm_wiring_sheet_name='herm chem synapse adjacency',  # herm gap jn synapse adjacency'
    #     name="2020_si_2_herm_chem_synapse",
    #     )

    # network = loader.load_network_file(
    #     worm_wiring_xlsx_file_path="networks/data/Cook_2019/SI 5 Connectome adjacency matrices, corrected July "
    #                                "2020.xlsx",
    #     worm_wiring_sheet_name='hermaphrodite chemical',
    #     name="2020_si_5_herm_chem_synapse",
    # )

    # network = loader.load_network_file(adj_file_path="networks/data/Uri_Alon_2002/example.txt",
    #                                    name="paper example")

    # network = loader.load_network_file(adj_file_path="networks/data/Uri_Alon_2002/coliInterNoAutoRegVec.txt",
    #                                    name='colinet1_noAuto')

    # network = loader.load_network_file(polarity_xlsx_file_path="networks/data/polarity_2020/s1_data.xlsx",
    #                                    polarity_sheet_name='5. Sign prediction',
    #                                    name="polarity 2020 SI 1")

    network = loader.load_network_file(durbin_file_path="networks/data/Durbin_1986/neurodata.txt",
                                       name='durbin network')

    # motif_search(network=network)
