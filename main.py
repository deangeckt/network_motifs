from networkx import DiGraph
from tqdm import tqdm

from motif_criteria import MotifCriteria
from network import Network
from network_randomizer import NetworkRandomizer
from subgraphs.induced_mfinder_enumeration import MFinderInduced
from subgraphs.none_induced_mfinder_enumeration import MFinderNoneInduced
from subgraphs.specific_subgraphs import SpecificSubGraphs
from subgraphs.sub_graphs import SubGraphs
from subgraphs.sub_graphs_utils import get_sub_id_name, get_sub_graph_from_id, SubGraphAlgoName, MotifName, \
    UniqueSubGraph
from utils.config import Config
from utils.simple_logger import Logger, LogLvl
import networkx as nx
import time

sub_graph_algorithms = {
    SubGraphAlgoName.specific: SpecificSubGraphs,
    SubGraphAlgoName.mfinder_induced: MFinderInduced,
    SubGraphAlgoName.mfinder_none_induced: MFinderNoneInduced
}


def sub_graph_search(network: DiGraph) -> dict:
    sub_graph_algo: SubGraphs = sub_graph_algorithms[algo](network)
    network_sub_graphs = sub_graph_algo.search_sub_graphs(k=k)

    if algo != SubGraphAlgoName.specific:
        total = 0
        for sub_id in network_sub_graphs:
            sub_graph = get_sub_graph_from_id(decimal=sub_id, k=k)
            amount = network_sub_graphs[sub_id]
            total += amount
            sub_name = get_sub_id_name(sub_id=sub_id, k=k)
            sub_name_log = f'\nSubgraph {sub_id}:'
            if sub_name:
                sub_name_log += f' - {sub_name}'
            logger.info(sub_name_log)

            logger.info(nx.adjacency_matrix(sub_graph).todense())
            logger.info(f'Appearances: {amount}')

        logger.info(f'\nMotif candidates: {len(network_sub_graphs)}')
        logger.info(f'Total appearances: {total}')

    return network_sub_graphs


def motif_search(file_path: str, name: str):
    logger.info(f'Motif search - {algo} algorithm, for {name} network:\n')
    with open(file_path, "r") as f:
        network = Network()
        network.load_from_txt(f.readlines())
        network.log_properties()

        start_time = time.time()
        network_sub_graphs = sub_graph_search(network.graph)
        end_time = time.time()
        logger.info(f'Sub graph search timer [S]: {round(end_time - start_time, 2)}')

        if not config.get_boolean_property('run_args', 'full_search'):
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
                sub_name_log += f' - {sub_name}'
            logger.info(sub_name_log)

            random_network_samples = [rand_network[sub_id] for rand_network in random_network_sub_graphs if
                                      sub_id in rand_network]
            MotifCriteria().is_motif(network_sub_graphs[sub_id], random_network_samples)


def debug_compare_to_uri():
    file_path = "../colinet-1.0/coliInterNoAutoRegVec.txt"
    with open(file_path, "r") as f:
        network = Network()
        network.load_from_txt(f.readlines())
        network.log_properties()

        specific_search = SpecificSubGraphs(network.graph, search=[MotifName.cascades])
        specific_search.search_sub_graphs(k=3)
        all_subs = specific_search.debug_hash_

    file_path = "../mfinder/mfinder1.1/e_coli_3_sub12.txt"
    paper_all_subs = set()
    with open(file_path, "r") as f:
        for line in f.readlines():
            splitted = line.split()
            x = int(splitted[0])
            y = int(splitted[1])
            z = int(splitted[2])

            sub_graph = ((x, y), (y, z))
            sub_graph = UniqueSubGraph(sub_graph)
            paper_all_subs.add(sub_graph)

    print()
    print(f'my: {len(all_subs)}')
    print(f'paper: {len(paper_all_subs)}')

    for sg in all_subs:
        if sg not in paper_all_subs:
            e1, e2 = sg.sub_graph
            print(network.graph.has_edge(*e1))
            print(network.graph.has_edge(*e2))
            print(sg)

    print()
    for sg in paper_all_subs:
        if sg not in all_subs:
            e1, e2 = sg.sub_graph
            print(network.graph.has_edge(*e1))
            print(network.graph.has_edge(*e2))
            print(sg)


if __name__ == "__main__":
    logger = Logger(LogLvl.info)
    config = Config()
    # debug_compare_to_uri()

    k = int(config.get_property('run_args', 'k'))
    algo = SubGraphAlgoName(config.get_property('run_args', 'sub_graph_algorithm'))

    # g = nx.DiGraph([(0, 4), (0, 5), (3,1), (3,2), (3,5), (3,4), (0,1)])
    # g = nx.DiGraph([(1, 0), (2, 0), (1, 2), (0, 0)])
    # sub_graph_search(g)

    # motif_search("networks/toy_network.txt", "toy")
    # motif_search("networks/paper_exmp_network.txt", "paper example")

    # restore paper result on e.coli (use ~100-200 rand networks)
    motif_search("networks/coliInterNoAutoRegVec.txt", "colinet1_noAuto")  # mfinder ~ k=3: 2.3 sec, k=4: 67 sec
