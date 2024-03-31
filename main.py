from networkx import DiGraph
from tqdm import tqdm

from motif_criteria import MotifCriteria
from network import Network
from network_randomizer import NetworkRandomizer
from subgraphs.mfinder_enum_induced import MFinderInduced
from subgraphs.mfinder_enum_none_induced import MFinderNoneInduced
from subgraphs.specific_subgraphs import SpecificSubGraphs
from subgraphs.sub_graphs_abc import SubGraphsABC
from subgraphs.sub_graphs_utils import get_sub_id_name, get_sub_graph_from_id, SubGraphAlgoName
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


def sub_graph_search(network: Network) -> dict:
    sub_graph_algo: SubGraphsABC = sub_graph_algorithms[algo](network.graph)

    logger.info(f'Sub graph search using k: {k}')
    start_time = time.time()
    network_sub_graphs = sub_graph_algo.search_sub_graphs(k=k)
    end_time = time.time()
    logger.info(f'Sub graph search timer [S]: {round(end_time - start_time, 2)}')

    network_sub_graphs_full = None
    if config.get_boolean_property('run_args', 'run_map_sub_graphs'):
        network_sub_graphs_full = sub_graph_algo.get_sub_graphs_fully_mapped()

    if algo != SubGraphAlgoName.specific:
        total = 0
        for sub_id in network_sub_graphs:
            sub_graph = get_sub_graph_from_id(decimal=sub_id, k=k)
            amount = network_sub_graphs[sub_id]
            total += amount
            sub_name = get_sub_id_name(sub_id=sub_id, k=k)
            sub_name_log = f'\nSubgraph {sub_id}:'
            if sub_name:
                sub_name_log += f' {sub_name}'

            logger.info(sub_name_log)
            logger.info(str(nx.adjacency_matrix(sub_graph).todense()))
            logger.info(f'Appearances: {amount}')

            if network_sub_graphs_full:
                logger.info(f'Appearances (edges): {network_sub_graphs_full[sub_id]}')
                sorted_nodes = network.sort_node_appearances_in_sub_graph(network_sub_graphs_full[sub_id])
                logger.info(f'Appearances Sorted by Nodes: {sorted_nodes}')

        logger.info(f'\nMotif candidates: {len(network_sub_graphs)}')
        logger.info(f'Total appearances: {total}')

    return network_sub_graphs


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

    network_sub_graphs = sub_graph_search(network)

    if not config.get_boolean_property('run_args', 'run_motif_criteria'):
        return

    randomizer = NetworkRandomizer(network.graph)
    random_network_amount = int(config.get_property('random', 'network_amount'))
    random_networks = randomizer.generate(amount=random_network_amount)

    logger.toggle(False)
    random_network_sub_graphs = [sub_graph_algorithms[algo](network).search_sub_graphs(k=k) for network in
                                 tqdm(random_networks)]
    logger.toggle(True)

    for sub_id in network_sub_graphs:
        sub_name = get_sub_id_name(sub_id=sub_id, k=k)
        sub_name_log = f'\nSubgraph {sub_id}:'
        if sub_name:
            sub_name_log += f' {sub_name}'
        logger.info(sub_name_log)

        random_network_samples = [rand_network[sub_id] for rand_network in random_network_sub_graphs if
                                  sub_id in rand_network]
        MotifCriteria().is_motif(network_sub_graphs[sub_id], random_network_samples)


if __name__ == "__main__":
    logger = Logger(LogLvl.info)
    config = Config()

    k = int(config.get_property('run_args', 'k'))
    algo = SubGraphAlgoName(config.get_property('run_args', 'sub_graph_algorithm'))

    # g = nx.DiGraph([(0, 4), (0, 5), (3,1), (3,2), (3,5), (3,4), (0,1)])
    # sub_graph_search(g)

    # motif_search("networks/toy_network.txt", "toy")
    # motif_search("networks/Uri_Alon_2002/example.txt", "paper example")

    # motif_search("networks/Uri_Alon_2002/coliInterNoAutoRegVec.txt", "colinet1_noAuto")

    motif_search("networks/Cook_2019/2020_si_2_herm_chem_synapse_adj.txt",
                 "2020_si2_herm_chem_synapse",
                 "networks/Cook_2019/2020_si_2_herm_neurons.txt")

    # motif_search("networks/Cook_2019/2020_si_2_herm_gap_adj.txt",
    #              "2020_si2_herm_gap",
    #              "networks/Cook_2019/2020_si_2_herm_neurons.txt")
