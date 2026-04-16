import pytest
import numpy as np
from src.fea.node import Node

def test_node_initialization():
    pos = np.zeros(6)
    node = Node(idt=1, position=pos)
    assert node.idt == 1
    np.testing.assert_array_equal(node.position, pos)
    np.testing.assert_array_equal(node.force, np.zeros(6))
    np.testing.assert_array_equal(node.displacements, np.zeros(6))

def test_node_distance():
    node1 = Node(idt=1, position=np.zeros(6))
    pos2 = np.zeros(6)
    pos2[0] = 3.0
    pos2[1] = 4.0
    node2 = Node(idt=2, position=pos2)
    assert node1.distance(node2) == 5.0

def test_node_dofs():
    node = Node(idt=1, position=np.zeros(6))
    # Com idt=1, dofs devem ser [1*6-1, 1*6-2, ..., 1*6-6] = [5, 4, 3, 2, 1, 0] ?
    # _calculate_dofs: np.array([self.idt * 6 - i for i in range(1, 7)])
    # Para idt=1: [6-1, 6-2, 6-3, 6-4, 6-5, 6-6] = [5, 4, 3, 2, 1, 0]
    np.testing.assert_array_equal(node.dofs, np.array([5, 4, 3, 2, 1, 0]))
    
    node2 = Node(idt=2, position=np.zeros(6))
    # Para idt=2: [12-1, 12-2, 12-3, 12-4, 12-5, 12-6] = [11, 10, 9, 8, 7, 6]
    np.testing.assert_array_equal(node2.dofs, np.array([11, 10, 9, 8, 7, 6]))

def test_node_displaced_position():
    pos = np.array([1.0, 2.0, 3.0, 0.0, 0.0, 0.0])
    node = Node(idt=1, position=pos)
    node.displacements = np.array([0.1, -0.2, 0.5, 0.0, 0.0, 0.0])
    
    expected_pos = np.array([1.1, 1.8, 3.5, 0.0, 0.0, 0.0])
    np.testing.assert_array_almost_equal(node.displaced_position(), expected_pos)
