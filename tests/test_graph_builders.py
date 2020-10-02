# Python packages
import unittest
import networkx as nx
import itertools as it
# Local modules
from gsc.utils import flatten
from gsc.graph_builders import linear_graph, make_crazy, from_MDS_code, \
    create_prime_graph, create_prime_power_graph


def process_graph_nodes_edges(graph, data=None):
    g_nodes = sorted(graph.nodes())
    if data is None:
        g_edges = sorted([tuple(sorted(edge)) for edge in graph.edges()])
    else:
        g_edges = sorted([tuple(sorted((u, v)) + [d]) for u, v, d
                          in graph.edges(data='weight')])
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

    def test_from_MDS_code(self):
        """ Tests building an AME graph state from an MDS code """
        prime, power = 5, 1
        A = [[1, 1, 1], [1, 2, 3], [1, 3, 4]]
        graph = from_MDS_code(A, prime, power)
        target_edges = [(0, 3, 1), (0, 4, 1), (0, 5, 1), (1, 3, 1), (1, 4, 2),
                        (1, 5, 3), (2, 3, 1), (2, 4, 3), (2, 5, 4)]
        self.assertEqual(type(graph), nx.Graph)
        self.assertEqual(list(graph.edges(data='weight')), target_edges)

    def test_create_prime_graph(self):
        nodes, prime = 10, 5
        w_edges = [(i, i+1 % nodes, i % prime) for i in range(nodes)]
        g = create_prime_graph(w_edges, prime)
        self.assertEqual(g.prime, prime)
        self.assertEqual(g.power, 1)
        self.assertEqual(g.dimension, prime)
        self.assertEqual(type(g), nx.Graph)
        self.assertEqual(list(g.edges(data='weight')), w_edges)

    def test_create_prime_power_graph(self):
        nodes, prime, power = 6, 7, 3
        all_nodes = [(i, j) for i in range(nodes) for j in range(power)]
        w_edges = [((0, 0), (1, 0), 1), ((0, 0), (1, 1), 3),
                   ((0, 0), (2, 0), 6), ((0, 0), (4, 2), 2),
                   ((0, 2), (5, 2), 4), ((1, 0), (3, 0), 5),
                   ((1, 0), (4, 2), 3), ((1, 1), (2, 2), 2),
                   ((2, 1), (4, 2), 2), ((4, 0), (5, 1), 5)]
        g = create_prime_power_graph(w_edges, prime, power)
        nodes, edges = process_graph_nodes_edges(g, data='weight')
        self.assertEqual(type(g), nx.Graph)
        self.assertEqual(g.prime, prime)
        self.assertEqual(g.power, power)
        self.assertEqual(edges, w_edges)
        self.assertEqual(nodes, all_nodes)


if __name__ == '__main__':
    unittest.main()
