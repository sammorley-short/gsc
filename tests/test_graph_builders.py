# Python packages
import pytest
import itertools as it
import networkx as nx
# Local modules
from gsc.utils import flatten
from gsc.graph_builders import (
    linear_graph,
    make_crazy,
    from_MDS_code,
    create_prime_graph,
    create_prime_power_graph,
)


@pytest.mark.parametrize('n_v', [8, 9, 10])
def test_linear_graph_builder(n_v):
    graph = linear_graph(n_v)
    nodes, edges = process_graph_nodes_edges(graph)
    expected_edges = sorted([((i, 0), (i + 1, 0)) for i in range(n_v - 1)])
    expected_nodes = sorted([(i, 0) for i in range(n_v)])
    assert nodes == expected_nodes
    assert edges == expected_edges


@pytest.mark.parametrize('length', [8, 9, 10])
@pytest.mark.parametrize('height', [1, 2, 3])
def test_crazy_linear_graph_builder(length, height):
    crazy_graph = make_crazy(linear_graph(length), height)
    nodes, edges = process_graph_nodes_edges(crazy_graph)
    encoded_nodes = {(i, 0): [((i, 0), j) for j in range(height)] for i in range(length)}
    encoded_edges = [((i, 0), (i + 1, 0)) for i in range(length - 1)]
    expected_edges = sorted(flatten([
        it.product(encoded_nodes[u], encoded_nodes[v])
        for u, v in encoded_edges
    ]))
    expected_nodes = sorted(flatten(encoded_nodes.values()))
    assert nodes == expected_nodes
    assert edges == expected_edges


@pytest.mark.parametrize('prime, power, A_matrix, expected_edges', [
    (
        5, 1,
        [
            [1, 1, 1],
            [1, 2, 3],
            [1, 3, 4]
        ],
        [
            (0, 3, 1), (0, 4, 1), (0, 5, 1), (1, 3, 1), (1, 4, 2),
            (1, 5, 3), (2, 3, 1), (2, 4, 3), (2, 5, 4)
        ]
    ),

])
def test_from_MDS_code(prime, power, A_matrix, expected_edges):
    """ Tests building an AME graph state from an MDS code """
    graph = from_MDS_code(A_matrix, prime, power)
    assert isinstance(graph, nx.Graph)
    assert list(graph.edges(data='weight')) == expected_edges


@pytest.mark.parametrize('n_v', [8, 9, 10])
@pytest.mark.parametrize('prime', [2, 3, 5])
def test_create_prime_graph(n_v, prime):
    w_edges = [(i, i + 1 % n_v, i % prime) for i in range(n_v)]
    graph = create_prime_graph(w_edges, prime)
    assert graph.prime == prime
    assert graph.power == 1
    assert graph.dimension == prime
    assert isinstance(graph, nx.Graph)
    assert list(graph.edges(data='weight')) == w_edges


@pytest.mark.parametrize('n_v, prime, power, weighted_edges', [
    (
        6, 7, 3, [
            ((0, 0), (1, 0), 1), ((0, 0), (1, 1), 3),
            ((0, 0), (2, 0), 6), ((0, 0), (4, 2), 2),
            ((0, 2), (5, 2), 4), ((1, 0), (3, 0), 5),
            ((1, 0), (4, 2), 3), ((1, 1), (2, 2), 2),
            ((2, 1), (4, 2), 2), ((4, 0), (5, 1), 5)
        ]
    ),
])
def test_create_prime_power_graph(n_v, prime, power, weighted_edges):
    all_nodes = [(i, j) for i in range(n_v) for j in range(power)]
    graph = create_prime_power_graph(weighted_edges, prime, power)
    nodes, edges = process_graph_nodes_edges(graph, data='weight')
    assert isinstance(graph, nx.Graph)
    assert graph.prime == prime
    assert graph.power == power
    assert edges == weighted_edges
    assert nodes == all_nodes


def process_graph_nodes_edges(graph, data=None):
    g_nodes = sorted(graph.nodes())
    if data is None:
        g_edges = sorted([tuple(sorted(edge)) for edge in graph.edges()])
    else:
        g_edges = sorted([
            tuple(sorted((u, v)) + [d])
            for u, v, d in graph.edges(data='weight')
        ])
    return g_nodes, g_edges
