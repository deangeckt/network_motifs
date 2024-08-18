import pickle
from utils.types import SearchResultBinaryFile, NetworkBinaryFile


def export_results(data: SearchResultBinaryFile):
    with open(data['args'].bin_file, 'wb') as f:
        pickle.dump(data, f)


def import_results(file_path: str) -> SearchResultBinaryFile:
    with open(file_path, 'rb') as f:
        return SearchResultBinaryFile(**pickle.load(f))


def export_network(file_path: str, data: NetworkBinaryFile):
    with open(file_path, 'wb') as f:
        pickle.dump(data, f)


def import_network(file_path: str) -> NetworkBinaryFile:
    with open(file_path, 'rb') as f:
        return NetworkBinaryFile(**pickle.load(f))
