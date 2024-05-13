import pandas as pd

from networks.loaders.network_loader_strategy import NetworkLoaderStrategy
from post_motif_analysis.polarity_counter import count_network_polarity_ratio
from utils.types import NetworkLoaderArgs


class NeuronalPolarityLoader(NetworkLoaderStrategy):
    def __init__(self, args: NetworkLoaderArgs):
        """
        xlsx files from the paper: Fenyves BG, Szilágyi GS, Vassy Z, Sőti C, Csermely P.
        Synaptic polarity and sign-balance prediction using gene expression data in the Caenorhabditis elegans chemical
        synapse neuronal connectome network
        """
        super().__init__(args)

        # polarity configuration
        self.src_col = 0
        self.tar_col = 3
        self.weight_col = 4
        self.edge_type_col = 5
        self.polarity_col = 16
        self.prim_nt_col = 1

        # polarity options [+, -, no pred, complex]
        self.filter_polarity: list[str] = args.filter_polarity
        # primary neurotransmitter options [GABA, Glu, ACh, 0] (0 is an int)
        self.filter_prim_nt: list[str] = args.filter_prim_nt

    def load(self, *args):
        xlsx_path, sheet_name = args
        xls = pd.ExcelFile(xlsx_path)
        df = xls.parse(sheet_name, header=None)

        # filter
        self.logger.info(f'Filtering Neurons with polarity: {self.filter_polarity}')
        df = df[df[self.polarity_col].isin(self.filter_polarity)]
        self.logger.info(f'Filtering Neurons with primary neurotransmitter: {self.filter_prim_nt}')
        df = df[df[self.prim_nt_col].isin(self.filter_prim_nt)]
        if not len(df):
            raise Exception('Filtering results with an empty data')

        src_neurons_names = df.iloc[:, self.src_col]
        tar_neurons_names = df.iloc[:, self.tar_col]
        edge_weights = df.iloc[:, self.weight_col]  # these are the amount of synapses
        polarity = df.iloc[:, self.polarity_col]

        self.network.use_polarity = True
        full_graph_polarity_ratio = count_network_polarity_ratio(polarity)
        self.logger.info(f'Polarity E/I ratio (before filtering): {round(full_graph_polarity_ratio, 3)}')

        self.network.neuron_names = list(set(src_neurons_names) | set(tar_neurons_names))
        # make it deterministic between runs
        self.network.neuron_names.sort()
        neurons_indices = {ss: i for i, ss in enumerate(self.network.neuron_names)}

        for v1, v2, w, p in zip(src_neurons_names, tar_neurons_names, edge_weights, polarity):
            self._load_synapse(neurons_indices[v1], neurons_indices[v2], w, p)

        polarities = [(self.network.graph.get_edge_data(s, t)['polarity']) for s, t in self.network.graph.edges]
        self.network.polarity_ratio = count_network_polarity_ratio(polarities)

