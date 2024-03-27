from enum import Enum
from typing import Optional

from singleton_decorator import singleton


class LogLvl(int, Enum):
    info = 0
    debug = 1


@singleton
class Logger:
    def __init__(self, lvl=LogLvl.info, file: Optional[str] = None):
        self.lvl = lvl
        self.on = True
        self.file = None
        if file is not None:
            self.file = open(file, "w")

    def change_file(self, file: str):
        if self.file is not None:
            self.file.close()
        self.file = open(file, "w")

    def toggle(self, on: bool):
        self.on = on

    def __print(self, msg: str):
        if not self.on:
            return
        if self.file is None:
            print(msg)
        else:
            self.file.write(msg+'\n')

    def debug(self, msg: str):
        if self.lvl < LogLvl.debug:
            return
        self.__print(msg)

    def info(self, msg: str):
        if self.lvl < LogLvl.info:
            return
        self.__print(msg)
