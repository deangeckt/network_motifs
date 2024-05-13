import pickle

from utils.types import BinaryFile


def export_results(data: BinaryFile):
    with open(data['args'].bin_file, 'wb') as f:
        pickle.dump(data, f)


def import_results(file_path: str) -> BinaryFile:
    with open(file_path, 'rb') as f:
        return BinaryFile(**pickle.load(f))
