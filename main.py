import random

from networkx import DiGraph
from tqdm import tqdm

from motif_criteria import MotifCriteria
from network import Network
from random_networks.markov_chain_switching import MarkovChainSwitching
from subgraphs.mfinder_enum_induced import MFinderInduced
from subgraphs.mfinder_enum_none_induced import MFinderNoneInduced
from subgraphs.specific_subgraphs import SpecificSubGraphs
from subgraphs.sub_graphs_abc import SubGraphsABC
from subgraphs.sub_graphs_utils import get_sub_id_name, get_sub_graph_from_id, SubGraphAlgoName, \
    generate_isomorphic_k_sub_graphs, get_number_of_disjoint_group_nodes
from utils.config import Config
from utils.simple_logger import Logger, LogLvl
import networkx as nx
import time
from typing import Optional

sub_graph_algorithms = {
    SubGraphAlgoName.specific: SpecificSubGraphs,
    SubGraphAlgoName.mfinder_induced: MFinderInduced,
    SubGraphAlgoName.mfinder_none_induced: MFinderNoneInduced
}


def custom_graph_search(graph: DiGraph):
    network = Network()
    network.load_graph(graph)
    sub_graph_search(network)


def log_motif_details(sub_id: int, network_sub_graphs: dict, network: Network, network_sub_graphs_full: dict):
    sub_name = get_sub_id_name(sub_id=sub_id, k=k)
    sub_name_log = f'\nSub graph {sub_id}:'
    if sub_name:
        sub_name_log += f' {sub_name}'
    sub_graph = get_sub_graph_from_id(decimal=sub_id, k=k)

    logger.info(sub_name_log)
    logger.info(str(nx.adjacency_matrix(sub_graph).todense()))

    if sub_id not in network_sub_graphs:
        return

    amount = network_sub_graphs[sub_id]
    logger.info(f'Appearances: {amount}')

    if config.get_boolean_property('run_args', 'log_all_motif_sub_graphs'):
        logger.info(f'Appearances (edges): {network_sub_graphs_full[sub_id]}')
        sorted_nodes = network.sort_node_appearances_in_sub_graph(network_sub_graphs_full[sub_id])
        logger.info(f'Appearances of Nodes in the sub-graph: {sorted_nodes}')


def sub_graph_search(network: Network) -> tuple[dict, dict]:
    sub_graph_algo: SubGraphsABC = sub_graph_algorithms[algo](network.graph, isomorphic_mapping)

    logger.info(f'Sub graph search using k: {k}')
    logger.info(f'Total num of different sub graphs size {k} is: {len(isomorphic_graphs)}')
    allow_self_loops = config.get_boolean_property('run_args', 'allow_self_loops')
    logger.info(f'Allow self loops: {allow_self_loops}')

    start_time = time.time()
    network_sub_graphs = sub_graph_algo.search_sub_graphs(k=k)
    end_time = time.time()
    logger.info(f'Sub graph search timer [Sec]: {round(end_time - start_time, 2)}')
    network_sub_graphs_full = sub_graph_algo.get_sub_graphs_fully_mapped()

    total_sub_graphs = sum(network_sub_graphs.values())
    logger.info(f'\nMotif candidates: {len(network_sub_graphs)}')
    logger.info(f'Total number of {k}-node sub graphs: {total_sub_graphs}')

    if algo != SubGraphAlgoName.specific:
        if not config.get_boolean_property('run_args', 'run_motif_criteria'):
            for sub_id in network_sub_graphs:
                log_motif_details(sub_id, network_sub_graphs, network, network_sub_graphs_full)
                logger.info(f'Uniq: {get_number_of_disjoint_group_nodes(network_sub_graphs_full[sub_id])}')

    return network_sub_graphs, network_sub_graphs_full


def motif_search(network: Network):
    logger.info(f'Algorithm: {algo}\n')

    if not config.get_boolean_property('run_args', 'run_sub_graph_search'):
        return

    network_sub_graphs, network_sub_graphs_full = sub_graph_search(network)

    if not config.get_boolean_property('run_args', 'run_motif_criteria'):
        return

    randomizer = MarkovChainSwitching(network.graph)
    random_network_amount = int(config.get_property('random', 'network_amount'))
    random_networks = randomizer.generate(amount=random_network_amount)

    logger.toggle(False)
    random_network_sub_graphs = []
    for rand_network in tqdm(random_networks):
        sub_graph_algo: SubGraphsABC = sub_graph_algorithms[algo](rand_network, isomorphic_mapping)
        random_network_sub_graphs.append(sub_graph_algo.search_sub_graphs(k=k))
    logger.toggle(True)

    motif_criteria = MotifCriteria()

    for sub_id in tqdm(isomorphic_graphs):
        log_motif_details(sub_id, network_sub_graphs, network, network_sub_graphs_full)

        random_network_samples = [rand_network.get(sub_id, 0) for rand_network in random_network_sub_graphs]
        motif_criteria.is_motif(n_real=network_sub_graphs.get(sub_id, 0),
                                random_network_samples=random_network_samples,
                                sub_graphs=network_sub_graphs_full[sub_id]
                                )


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
    logger = Logger(LogLvl.info, 'cook_100_wo_mutual_edges.txt')
    config = Config()

    k = int(config.get_property('run_args', 'k'))
    algo = SubGraphAlgoName(config.get_property('run_args', 'sub_graph_algorithm'))
    isomorphic_mapping, isomorphic_graphs = generate_isomorphic_k_sub_graphs(k=k)

    # custom_graph_search(nx.DiGraph([(1, 0), (2, 0), (1, 2)]))

    # network = load_network_file(adj_file_path="networks/toy_network.txt", name='toy')
    #
    # network = load_network_file(adj_file_path="networks/Uri_Alon_2002/example.txt", name="paper example")
    #
    # network = load_network_file(adj_file_path="networks/Uri_Alon_2002/coliInterNoAutoRegVec.txt",
    #                             name='colinet1_noAuto')
    #
    network = load_network_file(adj_file_path="networks/Cook_2019/2020_si_2_herm_chem_synapse_adj_5.txt",
                                neurons_file_path="networks/Cook_2019/2020_si_2_herm_neurons.txt",
                                name="2020_si2_herm_chem_synapse_5"
                                )

    # network = load_network_file(polarity_xlsx_file_path="networks/polarity_2020/s1_data.xlsx",
    #                             polarity_sheet_name='5. Sign prediction',
    #                             name="polarity 2020 SI 1")
    motif_search(network=network)
