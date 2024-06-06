import numpy as np
import pandas as pd

from networks.loaders.network_loader_strategy import NetworkLoaderStrategy
from utils.types import NetworkLoaderArgs

from pydantic import BaseModel


class XlsxParams(BaseModel):
    sheet_name: str
    row_start: int
    col_start: int
    name_start: int


class WormWiringLoader(NetworkLoaderStrategy):
    "Cook et al., 2019"

    def __init__(self, args: NetworkLoaderArgs):
        super().__init__(args)
        self.amount = 272

        self.filter_syn_type = args.filter_syn_type
        self.filter_sex_type = args.filter_sex_type

        self.logger.info(f'Filtering Neurons of: {self.filter_sex_type}')
        self.logger.info(f'Filtering Synapses of type: {self.filter_syn_type}')

        self.si2_herm_chem_params = XlsxParams(sheet_name='herm chem synapse adjacency',
                                               row_start=60,
                                               col_start=60,
                                               name_start=2)
        self.si2_herm_gap_params = XlsxParams(sheet_name='herm gap jn synapse adjacency',
                                              row_start=60,
                                              col_start=60,
                                              name_start=2)
        self.si2_male_chem_params = XlsxParams(sheet_name='male chem synapse adjacency',
                                               row_start=1,
                                               col_start=1,
                                               name_start=0)
        self.si2_male_gap_params = XlsxParams(sheet_name='male gap jn synapse adjacency',
                                              row_start=1,
                                              col_start=1,
                                              name_start=0)

        self.si5_herm_chem_params = XlsxParams(sheet_name='hermaphrodite chemical',
                                               row_start=23,
                                               col_start=53,
                                               name_start=2)
        self.si5_herm_gap_params = XlsxParams(sheet_name='hermaphrodite gap jn symmetric',
                                              row_start=60,
                                              col_start=60,
                                              name_start=2)
        self.si5_male_chem_params = XlsxParams(sheet_name='male chemical',
                                               row_start=23,
                                               col_start=53,
                                               name_start=2)
        self.si5_male_gap_params = XlsxParams(sheet_name='male gap jn symmetric',
                                              row_start=60,
                                              col_start=60,
                                              name_start=2)

    def __get_xlsx_params(self, xlsx_path: str, filter_syn_type: str) -> XlsxParams:
        if 'SI 2' in xlsx_path:
            if self.filter_sex_type == 'herm':
                return self.si2_herm_chem_params if filter_syn_type == 'chem' else self.si2_herm_gap_params
            elif self.filter_sex_type == 'male':
                return self.si2_male_chem_params if filter_syn_type == 'chem' else self.si2_male_gap_params
            else:
                raise Exception(f'Unsupported filter_syn_type: {filter_syn_type}')
        elif 'SI 5' in xlsx_path:
            if self.filter_sex_type == 'herm':
                return self.si5_herm_chem_params if filter_syn_type == 'chem' else self.si5_herm_gap_params
            elif self.filter_sex_type == 'male':
                return self.si5_male_chem_params if filter_syn_type == 'chem' else self.si5_male_gap_params
            else:
                raise Exception(f'Unsupported filter_syn_type: {filter_syn_type}')
        else:
            raise Exception(f'Unsupported xlsx file: {xlsx_path}')

    def __load_xlsx_sheet(self, xlsx_path: str, filter_syn_type: str):
        params = self.__get_xlsx_params(xlsx_path, filter_syn_type)

        xls = pd.ExcelFile(xlsx_path)
        df = xls.parse(params.sheet_name, header=None)

        adj_mat = np.array(
            df.iloc[params.row_start:params.row_start + self.amount, params.col_start:params.col_start + self.amount])
        neurons_col = list(df.iloc[params.row_start:params.row_start + self.amount, params.name_start])
        neurons_row = list(df.iloc[params.name_start, params.col_start:params.col_start + self.amount])

        if not neurons_row == neurons_col:
            raise Exception(f'invalid xlsx file or indices')

        return adj_mat, neurons_col

    def __adj_mat_to_graph(self, adj_mat: np.ndarray, syn_type: str):
        for i in range(len(adj_mat)):
            for j in range(len(adj_mat)):
                value = adj_mat[i, j]
                if np.isnan(value) or value == 0:
                    continue

                synapse = value if syn_type == 'chem' else 0
                gap = value if syn_type == 'gap' else 0
                self._append_edge(i, j, synapse=synapse, gap=gap)

    def load(self, *args):
        """
        https://wormwiring.org/pages/adjacency.html xlsx format
        Cook, (2019). Whole-animal connectomes of both Caenorhabditis elegans sexes.
        """
        xlsx_path = args[0]
        if self.filter_syn_type == 'all':
            chem_adj_mat, chem_neurons = self.__load_xlsx_sheet(xlsx_path, 'chem')
            gap_adj_mat, gap_neurons = self.__load_xlsx_sheet(xlsx_path, 'gap')

            if not gap_neurons == chem_neurons:
                raise Exception(f'invalid xlsx file or indices')

            self.network.neuron_names = gap_neurons
            self.__adj_mat_to_graph(chem_adj_mat, syn_type='chem')
            self.__adj_mat_to_graph(gap_adj_mat, syn_type='gap')
        else:
            adj_mat, neuron_names = self.__load_xlsx_sheet(xlsx_path, self.filter_syn_type)
            self.network.neuron_names = neuron_names
            self.__adj_mat_to_graph(adj_mat, syn_type=self.filter_syn_type)
