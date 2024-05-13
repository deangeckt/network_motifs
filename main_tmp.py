# TODO: mv to some analysis notebook

from motif_search_main import load_network_from_args
from utils.export_import import import_results
from utils.logs import log_motif_criteria_args, log_motif_results, log_sub_graph_args, log_randomizer_args, \
    log_polarity_motif_results
from utils.simple_logger import Logger
from utils.types import MotifCriteriaArgs

logger = Logger()

data = import_results('tst.bin')
args = data['args']
motifs = data['motifs']
polarity_motifs = data['polarity_motifs']
network = load_network_from_args(args)

log_motif_criteria_args(MotifCriteriaArgs(**vars(args)))
log_sub_graph_args(args)
log_randomizer_args(args)
log_motif_results(motifs, network)
log_polarity_motif_results(polarity_motifs)
print(3)