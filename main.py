from networkx import DiGraph
from tqdm import tqdm

from motif_criteria import MotifCriteria
from network import Network
from network_randomizer import NetworkRandomizer
from subgraphs.mfinder_enumeration import MFinder
from subgraphs.specific_subgraphs import SpecificSubGraphs
from subgraphs.sub_graphs_utils import get_sub_id_name, get_sub_graph_from_id
from utils.config import Config
from utils.simple_logger import Logger, LogLvl
import networkx as nx

logger = Logger(LogLvl.info)
config = Config()


def sub_graph_search(graph: DiGraph):
    pass


def motif_search(file_path: str, name: str):
    k = 3
    logger.info(f'Motif search for {name} network:\n')
    with open(file_path, "r") as f:
        network = Network()
        network.load_from_txt(f.readlines())
        network.log_properties()

        mfinder = MFinder(network.graph)
        network_sub_graphs = mfinder.search_sub_graphs(k=k)
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

        randomizer = NetworkRandomizer(network.graph)
        random_network_amount = int(config.get_property('random', 'network_amount'))
        random_networks = randomizer.generate(amount=random_network_amount)

        logger.toggle(False)
        random_network_sub_graphs = [MFinder(network).search_sub_graphs(k=k) for network in
                                     tqdm(random_networks)]
        logger.toggle(True)

        for sub_id in network_sub_graphs:
            logger.info(f"\nsubgraph - {sub_id}")
            random_network_samples = [rand_network[sub_id] for rand_network in random_network_sub_graphs if
                                      sub_id in rand_network]
            MotifCriteria().is_motif(network_sub_graphs[sub_id], random_network_samples)


def specific_motif_search(file_path: str, name: str):
    logger.info(f'Specific Motif search for {name} network:\n')
    with open(file_path, "r") as f:
        network = Network()
        network.load_from_txt(f.readlines())
        network.log_properties()
        # network.plot()

        randomizer = NetworkRandomizer(network.graph)
        random_network_amount = int(config.get_property('random', 'network_amount'))
        random_networks = randomizer.generate(amount=random_network_amount)

        network_sub_graphs = SpecificSubGraphs(network.graph).search_sub_graphs(k=3)
        logger.toggle(False)
        random_network_sub_graphs = [SpecificSubGraphs(network).search_sub_graphs(k=3) for network in
                                     tqdm(random_networks)]
        logger.toggle(True)

        for sub_graph in network_sub_graphs:
            logger.info(f"\nsubgraph - {sub_graph}")
            random_network_samples = [rand_network[sub_graph] for rand_network in random_network_sub_graphs]
            MotifCriteria().is_motif(network_sub_graphs[sub_graph], random_network_samples)

    logger.info('\n')


if __name__ == "__main__":
    # specific_motif_search("networks/toy_network.txt", "toy")
    motif_search("networks/paper_exmp_network.txt", "paper example")
    # specific_motif_search("networks/systems_biology_ex_network.txt", "uri alon course homework")

    # restore paper result on e.coli
    # specific_motif_search("../colinet-1.0/coliInterNoAutoRegVec.txt", "colinet1_noAuto")
