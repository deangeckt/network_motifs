from tqdm import tqdm

from motif_criteria import MotifCriteria
from network import Network
from network_randomizer import NetworkRandomizer
from subgraphs.mfinder_enum_induced import MFinderInduced
from subgraphs.mfinder_enum_none_induced import MFinderNoneInduced
from subgraphs.specific_subgraphs import SpecificSubGraphs
from subgraphs.sub_graphs_abc import SubGraphsABC
from subgraphs.sub_graphs_utils import get_sub_id_name, get_sub_graph_from_id, SubGraphAlgoName, \
    generate_isomorphic_k_sub_graphs
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
        logger.info(f'Appearances Sorted by Nodes: {sorted_nodes}')


def sub_graph_search(network: Network) -> tuple[dict, Optional[dict]]:
    sub_graph_algo: SubGraphsABC = sub_graph_algorithms[algo](network.graph, isomorphic_mapping)

    logger.info(f'Sub graph search using k: {k}')
    logger.info(f'Connected sub graphs - all options: {len(isomorphic_graphs)}')
    allow_self_loops = config.get_boolean_property('run_args', 'allow_self_loops')
    logger.info(f'Connected sub graphs - allow self loops: {allow_self_loops}')

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

    return network_sub_graphs, network_sub_graphs_full


def motif_search(file_path: str, name: str, neurons_file: Optional[str] = None):
    logger.info(f'Network name: {name}')
    logger.info(f'Algorithm: {algo}\n')

    network = Network()
    if neurons_file is not None:
        network.load_adj_file(file_path, is_synapse=True)
        network.load_neurons_file(neurons_file)
    else:
        network.load_adj_file(file_path, is_synapse=False)
    network.properties()

    if not config.get_boolean_property('run_args', 'run_sub_graph_search'):
        return

    network_sub_graphs, network_sub_graphs_full = sub_graph_search(network)

    if not config.get_boolean_property('run_args', 'run_motif_criteria'):
        return

    randomizer = NetworkRandomizer(network.graph)
    random_network_amount = int(config.get_property('random', 'network_amount'))
    random_networks = randomizer.generate(amount=random_network_amount)

    logger.toggle(False)
    sub_graph_algo: SubGraphsABC = sub_graph_algorithms[algo](network.graph, isomorphic_mapping)
    random_network_sub_graphs = [sub_graph_algo.search_sub_graphs(k=k) for network in tqdm(random_networks)]
    logger.toggle(True)

    motif_criteria = MotifCriteria()

    for sub_id in isomorphic_graphs:
        sub_name = get_sub_id_name(sub_id=sub_id, k=k)
        sub_name_log = f'\nSub graph {sub_id}:'
        if sub_name:
            sub_name_log += f' {sub_name}'
        logger.info(sub_name_log)

        random_network_samples = [rand_network.get(sub_id, 0) for rand_network in random_network_sub_graphs]
        motif_criteria.is_motif(n_real=network_sub_graphs.get(sub_id, 0),
                                random_network_samples=random_network_samples,
                                sub_graphs=network_sub_graphs_full[sub_id]
                                )

        log_motif_details(sub_id, network_sub_graphs, network, network_sub_graphs_full)


if __name__ == "__main__":
    logger = Logger(LogLvl.info)
    config = Config()

    k = int(config.get_property('run_args', 'k'))
    algo = SubGraphAlgoName(config.get_property('run_args', 'sub_graph_algorithm'))
    isomorphic_mapping, isomorphic_graphs = generate_isomorphic_k_sub_graphs(k=k)

    # g = nx.DiGraph([(0, 4), (0, 5), (3,1), (3,2), (3,5), (3,4), (0,1)])
    # sub_graph_search(g)

    # motif_search("networks/toy_network.txt", "toy")
    # motif_search("networks/Uri_Alon_2002/example.txt", "paper example")

    motif_search("networks/Uri_Alon_2002/coliInterNoAutoRegVec.txt", "colinet1_noAuto")

    # motif_search("networks/Cook_2019/2020_si_2_herm_chem_synapse_adj_5.txt",
    #              "2020_si2_herm_chem_synapse",
    #              "networks/Cook_2019/2020_si_2_herm_neurons.txt")

    # motif_search("networks/Cook_2019/2020_si_2_herm_gap_adj.txt",
    #              "2020_si2_herm_gap",
    #              "networks/Cook_2019/2020_si_2_herm_neurons.txt")
