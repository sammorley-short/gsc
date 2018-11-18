# Python packages
import ast
import numpy as np
from itertools import chain, combinations
from collections import defaultdict
from math import pi, cos, sin
from abp import GraphState
from abp.util import xyz
# Local modules


def canonical_edge_order(edges):
    return tuple(sorted(tuple(sorted(edge)) for edge in edges))


def circular_positions(nodes, r):
    """ Assigns circular coordinates to a set of crazy nodes """
    n = len(nodes)
    thetas = np.linspace(0, 2 * pi, n)[:-1]
    x_y_pos = [(nodes[0], (0, 0))]
    x_y_pos += [(node, (r*cos(theta), r*sin(theta)))
                for node, theta in zip(nodes[1:], thetas)]
    x_y_pos = [(node, xyz(*vector_add(node[0], shift)))
               for node, shift in x_y_pos]
    return x_y_pos


def to_GraphState(graph, r=0.2):
    """
    Converts a NetworkX graph into a GraphState.
    If graph is crazy, lays out crazy nodes radially.
    """
    # Defines qubits and coordinates for crazy or normal graph
    if hasattr(graph, 'crazy'):
        crazy_nodes = defaultdict(list)
        for node in graph.nodes():
            crazy_nodes[node[0]].append(node)
        nodes = []
        for crazy_node, sub_nodes in crazy_nodes.items():
            nodes += circular_positions(sorted(sub_nodes), r)
    else:
        nodes = [(node, xyz(*node)) for node in graph.nodes()]
    # Creates GraphState, adds qubits and acts CZs
    edges = graph.edges()
    gs = GraphState()
    for node, position in nodes:
        gs.add_qubit(node, position=position)
        gs.act_hadamard(node)
    gs.act_czs(*edges)
    return gs


def vector_scale(u, s):
    """ Scales a vector u by some scalar s """
    return (i * s for i in u)


def vector_add(u, v):
    """ Adds vectors u and v """
    return map(sum, zip(u, v))


def flatten(array, level=1):
    """ Flattens array to given level """
    for i in range(level):
        array = [item for sublist in array for item in sublist]
    return array


def powerset(s):
    return chain.from_iterable(combinations(s, r) for r in range(1, len(s) + 1))
