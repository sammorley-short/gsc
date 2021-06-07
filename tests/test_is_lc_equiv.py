# Python modules
import pytest
import networkx as nx
from random import randint
# Local modules
from gsc.graph_builders import random_connected_graph
from gsc.explore_lc_orbit import apply_qubit_LCs
from gsc.is_lc_equiv import get_adjacency_matrix, are_lc_equiv


@pytest.mark.parametrize('n_v', [4, 5, 6])
@pytest.mark.parametrize('n_lcs', [5, 10])
def test_is_lc_equiv(n_v, n_lcs):
    """
    Tests is_lc_equiv by creating graphs, applying LCs and checking for
    equivalence.
    """
    for _ in range(25):
        lc_nodes = [randint(0, n_v - 1) for _ in range(n_lcs)]
        graph_init = random_connected_graph(n_v)
        graph_fin = apply_qubit_LCs(graph_init, lc_nodes)
        lc_equiv, _ = are_lc_equiv(graph_init, graph_fin)
        assert lc_equiv


@pytest.mark.parametrize('n_v', [8, 9, 10])
def test_get_adjacency_matrix(n_v):
    """ Test get_adjacency_matrix against NetworkX version """
    for _ in range(1000):
        graph = random_connected_graph(n_v)
        nodes = sorted(graph.nodes())
        test_adj_mat, key = get_adjacency_matrix(graph)
        nx_adj_mat = nx.to_numpy_array(graph, nodelist=nodes, dtype=int)
        assert nodes == key
        assert test_adj_mat.tolist() == nx_adj_mat.tolist()
