from networkx import DiGraph
from tqdm import tqdm

from motif_criteria import MotifCriteria
from network import Network
from network_randomizer import NetworkRandomizer
from subgraphs.mfinder_enumeration import MFinder
from subgraphs.specific_subgraphs import SpecificSubGraphs
from subgraphs.sub_graphs import SubGraphs
from subgraphs.sub_graphs_utils import get_sub_id_name, get_sub_graph_from_id, SubGraphAlgo
from utils.config import Config
from utils.simple_logger import Logger, LogLvl
import networkx as nx


sub_graph_algorithms = {
    SubGraphAlgo.specific: SpecificSubGraphs,
    SubGraphAlgo.mfinder: MFinder,
}


def sub_graph_search(network: DiGraph) -> dict:
    sub_graph_algo: SubGraphs = sub_graph_algorithms[algo](network)
    network_sub_graphs = sub_graph_algo.search_sub_graphs(k=k)

    if algo != SubGraphAlgo.specific:
        total = 0
        for sub_id in network_sub_graphs:
            sub_graph = get_sub_graph_from_id(decimal=sub_id, k=k)
            amount = network_sub_graphs[sub_id]
            total += amount
            sub_name = get_sub_id_name(sub_id=sub_id, k=k)
            if sub_name:
                logger.info(f'\nSubgraph {sub_name} - {sub_id}:')
            else:
                logger.info(f'\nSubgraph {sub_id}:')
            logger.info(nx.adjacency_matrix(sub_graph).todense())
            logger.info(f'Appearances: {amount}')

        logger.info(f'\nMotif candidate amount: {len(network_sub_graphs)}')
        logger.info(f'Total appearances: {total}')

    return network_sub_graphs


def motif_search(file_path: str, name: str):
    logger.info(f'Motif search for {name} network:\n')
    with open(file_path, "r") as f:
        network = Network()
        network.load_from_txt(f.readlines())
        network.log_properties()

        network_sub_graphs = sub_graph_search(network.graph)

        randomizer = NetworkRandomizer(network.graph)
        random_network_amount = int(config.get_property('random', 'network_amount'))
        random_networks = randomizer.generate(amount=random_network_amount)

        logger.toggle(False)
        random_network_sub_graphs = [sub_graph_algorithms[algo](network).search_sub_graphs(k=k) for network in
                                     tqdm(random_networks)]
        logger.toggle(True)

        for sub_id in network_sub_graphs:
            logger.info(f"\nsubgraph - {sub_id}")
            random_network_samples = [rand_network[sub_id] for rand_network in random_network_sub_graphs if
                                      sub_id in rand_network]
            MotifCriteria().is_motif(network_sub_graphs[sub_id], random_network_samples)


if __name__ == "__main__":
    logger = Logger(LogLvl.info)
    config = Config()

    k = int(config.get_property('run_args', 'k'))
    algo = SubGraphAlgo(config.get_property('run_args', 'sub_graph_algorithm'))

    # g = nx.DiGraph([(0, 1), (1, 0), (0, 2), (1, 2), (2, 1), (2, 0), (0, 0)])
    # sub_graph_search(g)

    # motif_search("networks/toy_network.txt", "toy")
    motif_search("networks/paper_exmp_network.txt", "paper example")
    # motif_search("networks/systems_biology_ex_network.txt", "uri alon course homework")

    # restore paper result on e.coli
    # motif_search("../colinet-1.0/coliInterNoAutoRegVec.txt", "colinet1_noAuto")
