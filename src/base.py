import numpy as np

class Identifiable:
    def __init__(self, idt: int):
        self.idt = idt

    @property
    def idt(self):
        return self._idt

    @idt.setter
    def idt(self, value):
        if not isinstance(value, int):
            raise TypeError("idt must be an integer")
        elif value <= 0:
            raise ValueError("idt must be a positive integer")
        else:
            self._idt = value
