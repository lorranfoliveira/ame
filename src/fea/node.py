from ..base import Identifiable
import numpy as np
from numpy.typing import NDArray


class Node(Identifiable):
    def __init__(self, idt: int, position: NDArray[np.float64], force: NDArray[np.float64] = np.zeros(3)):
        super().__init__(idt)
        self.position = position
        self.force = force
        self.displacements = np.zeros(3)

    def distance(self, other: Node) -> float:
        return np.linalg.norm(self.position - other.position)

    def dofs(self) -> NDArray[np.int64]:
        return np.array([self.idt * 2 - 1, self.idt * 2])

    def displaced_position(self) -> NDArray[np.float64]:
        return self.position + self.displacements
