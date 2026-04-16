from ..base import Identifiable
import numpy as np


class Section(Identifiable):
    def __init__(self, idt):
        super().__init__(idt)


class RectangularSection(Section):
    def __init__(self, idt, b: float, h: float):
        super().__init__(idt)
        self.b = b
        self.h = h

        self._area = self._calculate_area()
        self._perimeter = self._calculate_perimeter()
        self._ix = self._calculate_ix()
        self._iy = self._calculate_iy()
        self._it = self._calculate_it()

    # -------------------------------- Data validation --------------------------------
    @property
    def b(self):
        return self._b

    @b.setter
    def b(self, value):
        if not isinstance(value, (int, float)):
            raise TypeError("b must be a number")
        elif value <= 0:
            raise ValueError("b must be a positive number")
        else:
            self._b = value

    @property
    def h(self):
        return self._h

    @h.setter
    def h(self, value):
        if not isinstance(value, (int, float)):
            raise TypeError("h must be a number")
        elif value <= 0:
            raise ValueError("h must be a positive number")
        else:
            self._h = value

    # -------------------------------- Methods --------------------------------
    def _calculate_area(self) -> float:
        return self.b * self.h

    def _calculate_perimeter(self) -> float:
        return 2 * (self.b + self.h)

    def _calculate_ix(self) -> float:
        return self.b * self.h ** 3 / 12

    def _calculate_iy(self) -> float:
        return self.h * self.b ** 3 / 12

    def _calculate_it(self, max_iter=10, tol: float = 1e-4) -> float:
        b = min(self.b, self.h)
        h = max(self.b, self.h)
        s = 0
        v1 = 0
        err = 2 * tol
        i = 1
        while err > tol or i < max_iter:
            s += 1 / (i ** 5) * np.tanh((i * np.pi * h) / (2 * b))
            v0 = v1
            v1 = (h * b ** 3 / 3) * (1 - 192 / (np.pi ** 5) * b / h * s)
            if i > 1:
                err = abs(v1 - v0) / v0
            i += 1
        return v1

    # -------------------------------- Properties --------------------------------
    @property
    def area(self) -> float:
        return self._area

    @property
    def perimeter(self) -> float:
        return self._perimeter

    @property
    def ix(self) -> float:
        return self._ix

    @property
    def iy(self) -> float:
        return self._iy

    @property
    def it(self) -> float:
        return self._it
