import ast

from singleton_decorator import singleton
import configparser


@singleton
class Config:
    def __init__(self):
        self.config = configparser.ConfigParser()
        self.config.read('utils/config.ini')

    def get_property(self, section: str, prop: str):
        return self.config[section][prop]

    def get_boolean_property(self, section: str, prop: str) -> bool:
        return self.config.getboolean(section, prop)

    def get_string_list(self, section: str, prop: str) -> list[str]:
        str_ = self.config[section][prop]
        return ast.literal_eval(str_)
