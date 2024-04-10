import random

from networkx import DiGraph
from tabulate import tabulate
from tqdm import tqdm

from motif_criteria import MotifCriteria
from network import Network
from random_networks.markov_chain_switching import MarkovChainSwitching
from subgraphs.mfinder_enum_induced import MFinderInduced
from subgraphs.mfinder_enum_none_induced import MFinderNoneInduced
from subgraphs.specific_subgraphs import SpecificSubGraphs
from subgraphs.sub_graphs_abc import SubGraphsABC
from subgraphs.sub_graphs_utils import generate_isomorphic_k_sub_graphs, create_base_motif
from utils.config import Config
from utils.simple_logger import Logger, LogLvl
import time
from typing import Optional

from utils.types import SubGraphAlgoName, Motif, MotifCriteriaResults

sub_graph_algorithms = {
    SubGraphAlgoName.specific: SpecificSubGraphs,
    SubGraphAlgoName.mfinder_induced: MFinderInduced,
    SubGraphAlgoName.mfinder_none_induced: MFinderNoneInduced
}


def custom_graph_search(graph: DiGraph):
    network = Network()
    network.load_graph(graph)
    sub_graph_search(network)
    # TODO: log?


def log_motif_details(sub_id: int, network_sub_graphs: dict, network: Network, network_sub_graphs_full: dict):
    if config.get_boolean_property('run_args', 'log_all_motif_sub_graphs'):
        logger.info(f'Appearances (edges): {network_sub_graphs_full[sub_id]}')
        sorted_nodes = network.sort_node_appearances_in_sub_graph(network_sub_graphs_full[sub_id])
        logger.info(f'Appearances of Nodes in the sub-graph: {sorted_nodes}')


def log_motif_results(motifs: dict[int, Motif]):
    table = []
    for motif_id in motifs:
        motif = motifs[motif_id]
        motif_criteria = motif.motif_criteria if motif.motif_criteria is not None else MotifCriteriaResults(
            n_real=motif.n_real,
        )

        table.append(
            [motif.id,
             motif.name,
             motif.adj_mat,
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
             ]
        )

    headers = ['ID', 'Name', 'Adj_mat', 'is-motif', 'N_real', 'N_rand', 'Std', 'Z_score', 'uniq',
               'is-significant', 'is-freq', 'is-uniq', 'is-anti-freq']

    col_align = tuple(['center'] * len(headers))
    float_fmt = tuple([".2f"] * len(headers))
    logger.info(tabulate(table, tablefmt="grid", headers=headers, colalign=col_align, floatfmt=float_fmt))


def sub_graph_search(network: Network) -> dict[int, Motif]:
    sub_graph_algo: SubGraphsABC = sub_graph_algorithms[algo](network.graph, isomorphic_mapping)

    logger.info(f'Sub Graph search using Algorithm: {algo}')
    logger.info(f'Sub Graph search using k: {k}')
    logger.info(f'Amount of isomorphic sub graphs of size k={k} is: {len(isomorphic_graphs)}')
    allow_self_loops = config.get_boolean_property('run_args', 'allow_self_loops')
    logger.info(f'Allow self loops: {allow_self_loops}')

    start_time = time.time()
    network_sub_graphs = sub_graph_algo.search_sub_graphs(k=k)
    end_time = time.time()
    logger.info(f'\nSub Graph search timer [Sec]: {round(end_time - start_time, 2)}')
    network_sub_graphs_full = sub_graph_algo.get_sub_graphs_fully_mapped()

    total_sub_graphs = sum(network_sub_graphs.values())
    logger.info(f'Motif candidates found: {len(network_sub_graphs)}')
    logger.info(f'Total number of {k}-node sub graphs found: {total_sub_graphs}')

    # TODO: handle specific algo and output
    # if algo != SubGraphAlgoName.specific:
    #     if not config.get_boolean_property('run_args', 'run_motif_criteria'):
    #         for sub_id in network_sub_graphs:
    #             log_motif_details(sub_id, network_sub_graphs, network, network_sub_graphs_full)
    #             logger.info(f'Uniq: {get_number_of_disjoint_group_nodes(network_sub_graphs_full[sub_id])}')

    motifs = {}
    for sub_id in network_sub_graphs:
        motif = create_base_motif(sub_id=sub_id, k=k)

        motif.n_real = network_sub_graphs[sub_id]
        motif.node_appearances = network.sort_node_appearances_in_sub_graph(network_sub_graphs_full[sub_id])
        motif.sub_graphs = network_sub_graphs_full[sub_id]

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
        random_network_samples = [rand_network.get(sub_id, 0) for rand_network in random_network_sub_graphs]

        motif_candidate: Motif = motif_candidates.get(sub_id, create_base_motif(sub_id=sub_id, k=k))
        motif_candidate.random_network_samples = random_network_samples
        motif_candidate.motif_criteria = motif_criteria.is_motif(motif_candidate)
        motif_candidates[sub_id] = motif_candidate

    log_motif_results(motif_candidates)


def load_network_file(
        name: str,
        adj_file_path: Optional[str] = None,
        neurons_file_path: Optional[str] = None,
        polarity_xlsx_file_path: Optional[str] = None,
        polarity_sheet_name: Optional[str] = None,
) -> Network:
    logger.info(f'Network name: {name}')
    network = Network()

    if adj_file_path and neurons_file_path:
        network.load_adj_neuronal_file(adj_file_path=adj_file_path, neurons_file_path=neurons_file_path)
    elif adj_file_path:
        network.load_adj_file(file_path=adj_file_path)
    elif polarity_xlsx_file_path and polarity_sheet_name:
        network.load_polarity_neuronal_file(xlsx_path=polarity_xlsx_file_path, sheet_name=polarity_sheet_name)
    else:
        logger.info('Error reading input file')
        exit(-1)

    network.properties()
    return network


if __name__ == "__main__":
    random.seed(42)
    logger = Logger(LogLvl.info)
    config = Config()

    k = int(config.get_property('run_args', 'k'))
    algo = SubGraphAlgoName(config.get_property('run_args', 'sub_graph_algorithm'))
    isomorphic_mapping, isomorphic_graphs = generate_isomorphic_k_sub_graphs(k=k)

    # custom_graph_search(nx.DiGraph([(1, 0), (2, 0), (1, 2)]))

    # network = load_network_file(adj_file_path="networks/toy_network.txt", name='toy')
    #
    network = load_network_file(adj_file_path="networks/Uri_Alon_2002/example.txt", name="paper example")
    #
    # network = load_network_file(adj_file_path="networks/Uri_Alon_2002/coliInterNoAutoRegVec.txt",
    #                             name='colinet1_noAuto')
    #
    # network = load_network_file(adj_file_path="networks/Cook_2019/2020_si_2_herm_chem_synapse_adj_5.txt",
    #                             neurons_file_path="networks/Cook_2019/2020_si_2_herm_neurons.txt",
    #                             name="2020_si2_herm_chem_synapse_5"
    #                             )

    # network = load_network_file(polarity_xlsx_file_path="networks/polarity_2020/s1_data.xlsx",
    #                             polarity_sheet_name='5. Sign prediction',
    #                             name="polarity 2020 SI 1")
    motif_search(network=network)
