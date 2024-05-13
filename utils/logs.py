from argparse import Namespace
from typing import Any

from networks.network import Network
from post_motif_analysis.polarity_counter import get_all_sub_graph_polarities
from utils.simple_logger import Logger
from utils.types import Motif, MotifCriteriaResults, MotifCriteriaArgs, SubGraphAlgoName, RandomGeneratorAlgoName
from tabulate import tabulate

logger = Logger()


def log_motifs_table(motifs: dict[Any, Motif]):
    if not motifs:
        return
    table = []
    for motif_id in motifs:
        motif = motifs[motif_id]
        motif_criteria_res = motif.motif_criteria if motif.motif_criteria is not None else MotifCriteriaResults(
            n_real=motif.n_real,
        )

        table.append([motif.id,
                      motif.name,
                      str(motif.adj_mat),
                      motif_criteria_res.is_motif,
                      motif_criteria_res.n_real,
                      motif_criteria_res.n_rand,
                      motif_criteria_res.std,
                      motif_criteria_res.z_score,
                      motif_criteria_res.uniq,
                      motif_criteria_res.is_statistically_significant,
                      motif_criteria_res.is_motif_frequent,
                      motif_criteria_res.is_uniq,
                      motif_criteria_res.is_anti_motif_frequent
                      ])

    headers = ['ID', 'Name', 'Adj_mat', 'is-motif', 'N_real', 'N_rand', 'Std', 'Z_score', 'uniq',
               'is-significant', 'is-freq', 'is-uniq', 'is-anti-freq']

    col_align = tuple(['center'] * len(headers))
    float_fmt = tuple([".2f"] * len(headers))
    logger.info(tabulate(table, tablefmt="grid", headers=headers, colalign=col_align, floatfmt=float_fmt))


def log_motif_results(motifs: dict[int, Motif], network: Network):
    motifs_list = motifs.values()
    non_zero_found = sum([1 for m in motifs_list if m.n_real > 0])
    total_sub_graphs = sum([m.n_real for m in motifs_list])

    logger.info(f'\nMotif candidates found: {non_zero_found}')
    logger.info(f'Total number of sub graphs found: {total_sub_graphs}')

    log_motifs_table(motifs)

    for motif_id in motifs:
        motif = motifs[motif_id]
        if motif.n_real == 0:
            continue
        logger.info(f'\nMotif Id: {motif_id}')

        # all_sub_graphs = get_all_sub_graph_polarities(motif.sub_graphs,
        #                                               network.graph) if network.use_polarity else motif.sub_graphs
        # logger.info(f'Appearances (all sub-graphs): {all_sub_graphs}')
        logger.info(f'Appearances of Nodes in the sub-graph: {motif.node_appearances}')
        logger.info(f'Role pattern: {motif.role_pattern}')
        for role in motif.node_roles:
            logger.info(f'Appearances of Nodes in role {role}: {motif.node_roles[role]}')

        for pol_freq in motif.polarity_frequencies:
            if pol_freq.frequency:
                logger.info(f'Polarity: {pol_freq.polarity} - frequency: {pol_freq.frequency}')


def log_polarity_motif_results(polarity_motifs: dict[int, dict[str, Motif]]):
    if polarity_motifs:
        logger.info('\nPolarity Motif results:')
        for polarity_motif in polarity_motifs.values():
            log_motifs_table(polarity_motif)


def log_motif_criteria_args(args: MotifCriteriaArgs):
    logger.info(f'\nMotif criteria:')
    logger.info(f'\talpha: {args.alpha}')
    logger.info(f'\tuse uniqueness: {args.use_uniq_criteria}')
    logger.info(f'\tuniqueness threshold: {args.uniqueness_threshold}')
    logger.info(f'\tfrequency threshold: {args.frequency_threshold}')


def log_sub_graph_args(args: Namespace):
    logger.info(f'\nSub Graph search using Algorithm: {SubGraphAlgoName(args.sub_graph_algorithm)}')
    logger.info(f'Sub Graph search using k: {args.k}')
    logger.info(f'Allow self loops: {args.allow_self_loops}')
    logger.info(f'Doing full Motif-Search: {args.run_motif_criteria}')


def log_randomizer_args(args: Namespace):
    random_generator_algo_choice = RandomGeneratorAlgoName(args.randomizer)
    logger.info(f'\nRandomizer: using {random_generator_algo_choice} algorithm')
    logger.info(f'Randomizer: generating {args.network_amount} random networks')
    if random_generator_algo_choice == RandomGeneratorAlgoName.markov_chain_switching:
        logger.info(f'Markov chain switch factor: {args.switch_factor}')

