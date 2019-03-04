# Python modules
import sys
import unittest
import itertools as it
import networkx as nx
from random import randint
# Local modules
from gsc.utils import canonical_edge_order
from gsc.graph_builders import random_connected_graph
from gsc.explore_lc_orbit import apply_qubit_LCs
from gsc.is_lc_equiv import get_adjacency_matrix, are_lc_equiv


class TestIsLCEquiv(unittest.TestCase):

    def setUp(self):
        pass

    def test_is_lc_equiv(self):
        """
            Tests is_lc_equiv by creating graphs, applying LCs and checking for
            equivalence.
        """
        n, lcs = 6, 10
        for _ in range(100):
            lc_nodes = [randint(0, n - 1) for _ in range(lcs)]
            graph_init = random_connected_graph(n)
            graph_fin = apply_qubit_LCs(graph_init, lc_nodes)
            lc_equiv, lc_ops = are_lc_equiv(graph_init, graph_fin)
            self.assertTrue(lc_equiv)

    def test_get_adjacency_matrix(self):
        """ Test get_adjacency_matrix against NetworkX version """
        n = 10
        for _ in range(1000):
            graph = random_connected_graph(n)
            nodes = sorted(graph.nodes())
            test_adj_mat, key = get_adjacency_matrix(graph)
            nx_adj_mat = nx.to_numpy_array(graph, nodelist=nodes, dtype=int)
            self.assertEqual(nodes, key)
            self.assertEqual(test_adj_mat.tolist(), nx_adj_mat.tolist())


if __name__ == '__main__':
    unittest.main()
