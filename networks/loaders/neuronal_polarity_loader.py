import pandas as pd

from networks.loaders.network_loader_strategy import NetworkLoaderStrategy
from utils.polarity_counter import count_network_polarity_ratio
from utils.types import NetworkLoaderArgs


class NeuronalPolarityLoader(NetworkLoaderStrategy):
    def __init__(self, args: NetworkLoaderArgs):
        """
        xlsx files from the paper: Fenyves BG, Szilágyi GS, Vassy Z, Sőti C, Csermely P.
        "Synaptic polarity and sign-balance prediction using gene expression data in the Caenorhabditis elegans chemical
        synapse neuronal connectome network"
        """
        super().__init__(args)

        # polarity configuration
        self.src_col = 0
        self.tar_col = 3
        self.weight_col = 4
        self.polarity_col = 16

        self.prim_nt_col = 1
        self.sec_nt_col = 2

        self.nt_ex_cols = {
            'Glu': 6,
            'ACh': 8,
            'GABA': 10
        }

        # polarity options [+, -, no pred, complex]
        self.filter_polarity: list[str] = args.filter_polarity
        # primary neurotransmitter options [GABA, Glu, ACh, 0] (0 is an int)
        self.filter_prim_nt: list[str] = args.filter_prim_nt

        self.sheet_name = '5. Sign prediction'

    @staticmethod
    def __get_nt_polarity_attr_aux(row: pd.DataFrame, col: int, nt: str, polarity: str):
        active = row[col]
        if active:
            return f'{nt}({polarity})'
        return ''

    def __get_nt_polarity_attr(self, row: pd.DataFrame, nt: str, polarity: str) -> str:
        if nt not in self.nt_ex_cols:
            return ''

        ex_attr = self.__get_nt_polarity_attr_aux(row, self.nt_ex_cols[nt], nt, '+')
        inh_attr = self.__get_nt_polarity_attr_aux(row, self.nt_ex_cols[nt] + 1, nt, '-')

        if polarity == 'complex':
            return f'{ex_attr} {inh_attr}'
        elif polarity == '+':
            return ex_attr
        else:
            return inh_attr

    def __get_full_polarity_attr(self, row: pd.DataFrame, nt1: str, nt2: str, polarity: str) -> str:
        nt1_attr = self.__get_nt_polarity_attr(row, nt1, polarity).strip()
        nt2_attr = self.__get_nt_polarity_attr(row, nt2, polarity).strip()
        return f'{nt1_attr} {nt2_attr}'.strip()

    def load(self, *args):
        xlsx_path = args[0]
        xls = pd.ExcelFile(xlsx_path)
        df = xls.parse(self.sheet_name, header=None)

        self.logger.info(f'Filtering Neurons with polarity: {self.filter_polarity}')
        df = df[df[self.polarity_col].isin(self.filter_polarity)]
        self.logger.info(f'Filtering Neurons with primary neurotransmitter: {self.filter_prim_nt}')
        df = df[df[self.prim_nt_col].isin(self.filter_prim_nt)]
        if not len(df):
            raise Exception('Filtering results with an empty data')

        src_neurons_names = df.iloc[:, self.src_col]
        tar_neurons_names = df.iloc[:, self.tar_col]
        edge_weights = df.iloc[:, self.weight_col]
        polarity = df.iloc[:, self.polarity_col]
        prim_nt = df.iloc[:, self.prim_nt_col]
        sec_nt = df.iloc[:, self.sec_nt_col]

        self.network.use_polarity = True
        full_graph_polarity_ratio = count_network_polarity_ratio(polarity)
        self.logger.info(f'Polarity ratios (before filtering): {full_graph_polarity_ratio}')

        self.network.neuron_names = list(set(src_neurons_names) | set(tar_neurons_names))
        # make it deterministic between runs
        self.network.neuron_names.sort()
        neurons_indices = {ss: i for i, ss in enumerate(self.network.neuron_names)}

        for idx, (v1, v2, w, p, nt1, nt2) in enumerate(zip(src_neurons_names,
                                                           tar_neurons_names,
                                                           edge_weights, polarity, prim_nt, sec_nt)):
            n1 = neurons_indices[v1]
            n2 = neurons_indices[v2]
            self._append_edge(n1, n2, synapse=w, gap=0, polarity=p)

            if self.network.graph.has_edge(n1, n2):
                self.network.graph[n1][n2]['full_polarity'] = self.__get_full_polarity_attr(df.iloc[idx], nt1, nt2, p)

        polarities = [(self.network.graph.get_edge_data(s, t)['polarity']) for s, t in self.network.graph.edges]
        self.network.polarity_ratio = count_network_polarity_ratio(polarities)
        self.network.polarity_options = self.filter_polarity
