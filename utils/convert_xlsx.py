import numpy as np
import pandas as pd


def convert_xlsx(xlsx_path: str, sheet_name: str, output_name: str):
    """
    2020 paper format....
    convert to the simple txt format: (v1, v2, w) per line.
    in this case the w isn't the weight but the number of synapses
    """
    xls = pd.ExcelFile(xlsx_path)
    df = xls.parse(sheet_name, header=None)

    start_idx = 60
    end_idx = 332
    names_idx = 2

    adj_mat = df.iloc[start_idx:end_idx, start_idx:end_idx]
    neurons_col = list(df.iloc[start_idx:end_idx, names_idx])
    neurons_row = list(df.iloc[names_idx, start_idx:end_idx])
    if not neurons_row == neurons_col:
        print(f'invalid xlsx file or indices')
        return

    neurons = neurons_col

    adj_output_file_name = f'{output_name}_adj.txt'
    with open(adj_output_file_name, "w") as f:
        amount_of_synapses = 0
        amount_of_edges = 0
        adj_mat = np.array(adj_mat)
        for i in range(len(adj_mat)):
            for j in range(len(adj_mat)):
                synapses = adj_mat[i, j]
                if not np.isnan(synapses):
                    amount_of_synapses += synapses
                    amount_of_edges += 1
                    f.write(f'{i} {j} {synapses}\n')

    neuron_output_file_name = f'{output_name}_neurons.txt'
    with open(neuron_output_file_name, "w") as f:
        [f.write(f'{neuron}\n') for neuron in neurons]

    print(f'total number of neurons: {len(neurons)}')
    print(f'total number of synapses: {amount_of_synapses}')
    print(f'total number of network edges: {amount_of_edges}')

    print(f'saved adj matrix file to: {adj_output_file_name}')
    print(f'saved neurons names file to: {neuron_output_file_name}')


if __name__ == "__main__":
    # convert_xlsx('nnn.xlsx', 'herm chem synapse adjacency', '2020_herm_chem_synapse')
    # convert_xlsx('nnn.xlsx', 'herm gap jn synapse adjacency', '2020_herm_gap')
    pass
