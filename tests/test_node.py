import pytest
import numpy as np
from src.fea.node import Node

def test_node_initialization():
    pos = np.array([0.0, 0.0, 0.0])
    node = Node(idt=1, position=pos)
    assert node.idt == 1
    np.testing.assert_array_equal(node.position, pos)
    np.testing.assert_array_equal(node.force, np.zeros(3))
    np.testing.assert_array_equal(node.displacements, np.zeros(3))

def test_node_distance():
    node1 = Node(idt=1, position=np.array([0.0, 0.0, 0.0]))
    node2 = Node(idt=2, position=np.array([3.0, 4.0, 0.0]))
    assert node1.distance(node2) == 5.0

def test_node_dofs():
    node = Node(idt=1, position=np.array([0.0, 0.0, 0.0]))
    # Com idt=1, dofs devem ser [1*2-1, 1*2] = [1, 2]
    np.testing.assert_array_equal(node.dofs(), np.array([1, 2]))
    
    node2 = Node(idt=2, position=np.array([1.0, 1.0, 1.0]))
    # Com idt=2, dofs devem ser [2*2-1, 2*2] = [3, 4]
    np.testing.assert_array_equal(node2.dofs(), np.array([3, 4]))

def test_node_displaced_position():
    pos = np.array([1.0, 2.0, 3.0])
    node = Node(idt=1, position=pos)
    node.displacements = np.array([0.1, -0.2, 0.5])
    
    expected_pos = np.array([1.1, 1.8, 3.5])
    np.testing.assert_array_almost_equal(node.displaced_position(), expected_pos)
