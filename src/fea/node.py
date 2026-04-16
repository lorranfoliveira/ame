from ..base import Identifiable
import numpy as np
from numpy.typing import NDArray


class Node(Identifiable):
    def __init__(self, idt: int, position: NDArray[np.float64], force: NDArray[np.float64] = np.zeros(6)):
        super().__init__(idt)
        self.position = position
        self.force = force
        self.displacements = np.zeros(6)

        self._dofs = self._calculate_dofs()

    @property
    def position(self) -> NDArray[np.float64]:
        return self._position

    @position.setter
    def position(self, value: NDArray[np.float64]):
        if not isinstance(value, np.ndarray):
            raise TypeError("Position must be a numpy array")
        elif value.shape != (6,):
            raise ValueError("Position must have 6 elements")
        else:
            self._position = value

    @property
    def force(self) -> NDArray[np.float64]:
        return self._force

    @force.setter
    def force(self, value: NDArray[np.float64]):
        if not isinstance(value, np.ndarray):
            raise TypeError("Force must be a numpy array")
        elif value.shape != (6,):
            raise ValueError("Force must have 6 elements")
        else:
            self._force = value

    @property
    def displacements(self) -> NDArray[np.float64]:
        return self._displacements

    @displacements.setter
    def displacements(self, value: NDArray[np.float64]):
        if not isinstance(value, np.ndarray):
            raise TypeError("Displacements must be a numpy array")
        elif value.shape != (6,):
            raise ValueError("Displacements must have 6 elements")
        else:
            self._displacements = value

    @property
    def dofs(self) -> NDArray[np.int64]:
        return self._dofs

    # -------------------------------- Methods --------------------------------
    def _calculate_dofs(self) -> NDArray[np.int64]:
        return np.array([self.idt * 6 - i for i in range(1, 7)])

    def distance(self, other: Node) -> float:
        return np.linalg.norm(self.position - other.position)

    def displaced_position(self, scale: float = 1) -> NDArray[np.float64]:
        return self.position + scale * self.displacements
