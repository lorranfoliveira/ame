from ..base import Identifiable


class Material(Identifiable):
    def __init__(self, idt: int, young: float, poisson: float, density: float):
        super().__init__(idt)
        self.young = young
        self.poisson = poisson
        self.density = density

    @property
    def young(self):
        return self._young

    @young.setter
    def young(self, value):
        if not isinstance(value, (int, float)):
            raise TypeError("Young's modulus must be a number")
        elif value <= 0:
            raise ValueError("Young's modulus must be a positive number")
        else:
            self._young = value

    @property
    def poisson(self):
        return self._poisson

    @poisson.setter
    def poisson(self, value):
        if not isinstance(value, (int, float)):
            raise TypeError("Poisson's ratio must be a number")
        elif value < 0 or value > 0.5:
            raise ValueError("Poisson's ratio must be between 0 and 0.5")
        else:
            self._poisson = value

    @property
    def density(self):
        return self._density

    @density.setter
    def density(self, value):
        if not isinstance(value, (int, float)):
            raise TypeError("Density must be a number")
        elif value <= 0:
            raise ValueError("Density must be a positive number")
        else:
            self._density = value

    @property
    def shear_modulus(self):
        return self.young / (2 * (1 + self.poisson))

    @property
    def bulk_modulus(self):
        return self.young / (3 * (1 - 2 * self.poisson))

    def to_dict(self) -> dict:
        return {
            "idt": self.idt,
            "young": self.young,
            "poisson": self.poisson,
            "density": self.density,
            "shear_modulus": self.shear_modulus,
            "bulk_modulus": self.bulk_modulus
        }
