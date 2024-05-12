from typing import Any

from networks.network import Network
from post_motif_analysis.polarity_counter import get_all_sub_graph_polarities
from utils.simple_logger import Logger
from utils.types import Motif, MotifCriteriaResults
from tabulate import tabulate

logger = Logger()


def log_motifs_table(motifs: dict[Any, Motif]):
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
    log_motifs_table(motifs)

    for motif_id in motifs:
        motif = motifs[motif_id]
        if motif.n_real == 0:
            continue
        logger.info(f'\nMotif Id: {motif_id}')

        all_sub_graphs = get_all_sub_graph_polarities(motif.sub_graphs,
                                                      network.graph) if network.use_polarity else motif.sub_graphs

        # logger.info(f'Appearances (all sub-graphs): {all_sub_graphs}')
        logger.info(f'Appearances of Nodes in the sub-graph: {motif.node_appearances}')
        logger.info(f'Role pattern: {motif.role_pattern}')
        for role in motif.node_roles:
            logger.info(f'Appearances of Nodes in role {role}: {motif.node_roles[role]}')

        for pol_freq in motif.polarity_frequencies:
            if pol_freq.frequency:
                logger.info(f'Polarity: {pol_freq.polarity} - frequency: {pol_freq.frequency}')
