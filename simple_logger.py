from enum import Enum
from singleton_decorator import singleton


class LogLvl(int, Enum):
    info = 0
    debug = 1


@singleton
class Logger:
    def __init__(self, lvl: LogLvl):
        self.lvl = lvl

    @staticmethod
    def __print(msg: str):
        print(msg)

    def debug(self, msg: str):
        if self.lvl < LogLvl.debug:
            return
        self.__print(msg)

    def info(self, msg: str):
        if self.lvl < LogLvl.info:
            return
        self.__print(msg)
