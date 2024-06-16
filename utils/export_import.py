import pickle

from utils.types import SearchResultBinaryFile


def export_results(data: SearchResultBinaryFile):
    with open(data['args'].bin_file, 'wb') as f:
        pickle.dump(data, f)


def import_results(file_path: str) -> SearchResultBinaryFile:
    with open(file_path, 'rb') as f:
        return SearchResultBinaryFile(**pickle.load(f))
