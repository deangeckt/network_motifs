from tqdm import tqdm

from motif_criteria import MotifCriteria
from network import Network
from network_randomizer import NetworkRandomizer
from subgraphs.mfinder_enumeration import MFinder
from subgraphs.specific_subgraphs import SpecificSubGraphs
from subgraphs.sub_graphs import get_sub_id_name
from utils.config import Config
from utils.simple_logger import Logger, LogLvl
import networkx as nx

logger = Logger(LogLvl.debug)
config = Config()


def motif_search(file_path: str, name: str):
    logger.info(f'Motif search for {name} network:\n')
    with open(file_path, "r") as f:
        network = Network()
        network.load_from_txt(f.readlines())
        network.log_properties()

        mfinder = MFinder(network.graph)
        network_sub_graphs = mfinder.search_sub_graphs(k=3)
        total = 0
        for sub_id in network_sub_graphs:
            sub_graph = mfinder.get_sub_graph_from_id(sub_id)
            amount = network_sub_graphs[sub_id]
            total += amount
            sub_name = get_sub_id_name(sub_id)
            if sub_name:
                print(f'\nSubgraph {sub_name} - {sub_id}:')
            else:
                print(f'\nSubgraph {sub_id}:')
            print(nx.adjacency_matrix(sub_graph).todense())
            print(f'Appearances: {amount}')

        print('\nMotif candidate amount:', len(network_sub_graphs))
        print('Total appearances:', total)


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

        network_sub_graphs = MFinder(network.graph).search_sub_graphs(k=3)
        logger.toggle(False)
        random_network_sub_graphs = [MFinder(network).search_sub_graphs(k=3) for network in
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
