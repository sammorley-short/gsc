# Python packages
import sys
import unittest
import networkx as nx
import itertools as it
# Local modules
sys.path.append('..')
from utils import flatten
from graph_builders import linear_graph, make_crazy


def process_graph_nodes_edges(graph):
    g_nodes = sorted(graph.nodes())
    g_edges = sorted([tuple(sorted(edge)) for edge in graph.edges()])
    return g_nodes, g_edges


class TestGraphBuilders(unittest.TestCase):

    def setUp(self):
        pass

    def test_linear_graph_builder(self):
        """ Tests a 10-node linear graph correctly built """
        l = 10
        g = linear_graph(l)
        g_nodes, g_edges = process_graph_nodes_edges(g)
        edges = sorted([((i, 0), (i+1, 0)) for i in range(l-1)])
        nodes = sorted([(i, 0) for i in range(l)])
        self.assertEqual(g_nodes, nodes)
        self.assertEqual(g_edges, edges)

    def test_crazy_linear_graph_builder(self):
        """ Tests a 10-node crazy linear graph correctly built """
        l, n = 5, 3
        g = linear_graph(l)
        cg = make_crazy(g, n)
        cg_nodes, cg_edges = process_graph_nodes_edges(cg)
        nodes = {(i, 0): [((i, 0), j) for j in range(n)] for i in range(l)}
        edges = [((i, 0), (i+1, 0)) for i in range(l-1)]
        edges = sorted(flatten([it.product(nodes[u], nodes[v])
                                for u, v in edges]))
        nodes = sorted(flatten(nodes.values()))
        self.assertEqual(cg_nodes, nodes)
        self.assertEqual(cg_edges, edges)


if __name__ == '__main__':
    unittest.main()
