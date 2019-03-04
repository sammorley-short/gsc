# Python packages
import sys
import random
import unittest
import networkx as nx
from abp import GraphState
from pprint import pprint
# Local modules
from gsc.utils import canonical_edge_order
from gsc.is_lc_equiv import are_lc_equiv
from gsc.explore_lc_orbit import qubit_LC, explore_lc_orbit
from gsc.graph_builders import create_prime_power_graph


def to_GraphState(graph):
    nodes = graph.nodes()
    edges = graph.edges()
    gs = GraphState()
    for node in nodes:
        gs.add_qubit(node)
        gs.act_hadamard(node)
    gs.act_czs(*edges)
    return gs


def gen_random_connected_graph(n, p=0.333):
    """ Generates a Erdos-Renyi G_n,p random graph """
    g = nx.fast_gnp_random_graph(n, p)
    while not nx.is_connected(g):
        g = nx.fast_gnp_random_graph(n, p)
    return g


class TestExploreClass(unittest.TestCase):

    def setUp(self):
        pass

    def test_qubit_LC(self):
        """ Tests local complementation works against abp version """
        for _ in range(100):
            # Creates a random NetworkX graph and it's equivalent GraphState
            g = gen_random_connected_graph(15)
            gs = to_GraphState(g)
            # Gets the original edges
            g_edges = canonical_edge_order(g.edges())
            gs_edges = canonical_edge_order(gs.edgelist())
            # Randomly picks a node for local complementation
            lc_node = random.choice(list(g.nodes()))
            # Performs local complementation on both graphs
            lc_g = qubit_LC(g, lc_node, 1)
            lc_gs = to_GraphState(g)
            lc_gs.local_complementation(lc_node)
            # Checks that their edgelists are equal
            lc_g_edges = canonical_edge_order(lc_g.edges())
            lc_gs_edges = canonical_edge_order(lc_gs.edgelist())
            self.assertEqual(lc_g_edges, lc_gs_edges)
            # Performs a second local complementation on same node
            lc_lc_g = qubit_LC(lc_g, lc_node, 1)
            lc_lc_gs = to_GraphState(lc_g)
            lc_lc_gs.local_complementation(lc_node)
            # Checks edgelists are equal and also equal to the originals
            lc_lc_g_edges = canonical_edge_order(lc_lc_g.edges())
            lc_lc_gs_edges = canonical_edge_order(lc_lc_gs.edgelist())
            self.assertEqual(lc_lc_g_edges, lc_lc_gs_edges)
            self.assertEqual(lc_lc_g_edges, g_edges)
            self.assertEqual(lc_lc_gs_edges, gs_edges)

    def test_explore_lc_orbit(self):
        """ Generates class graphs for random graphs """
        for _ in range(10):
            # Creates a random NetworkX graph and it's equivalent GraphState
            g = gen_random_connected_graph(7)
            class_graph = explore_lc_orbit(g, verbose=False)
            orbit_hashes = set([graph['hash'] for graph
                                in class_graph.node.values()])
            self.assertEqual(len(class_graph.node), len(orbit_hashes))
            graphs = [graph['nx_graph'] for graph in class_graph.node.values()]
            for i in range(10):
                graph_a = random.choice(graphs)
                graph_b = random.choice(graphs)
                lc_equiv, lc_ops = are_lc_equiv(graph_a, graph_b)
                self.assertTrue(lc_equiv)

    def test_ququart_pair(self):
        """ Tests working for ququart entangled pair LC classes """
        # Tests first equivalence class
        edges = [((0, 0), (1, 0), 1), ((0, 1), (1, 1), 1)]
        graph = create_prime_power_graph(edges, 2, 2)
        class_graph = explore_lc_orbit(graph, verbose=False)
        register = set(tuple(map(tuple, attrs['edges']))
                       for node, attrs in class_graph.node.iteritems())
        target = \
            [(((0, 1), (1, 0), 1), ((0, 0), (1, 1), 1)),
             (((0, 1), (1, 0), 1), ((0, 1), (1, 1), 1), ((0, 0), (1, 1), 1)),
             (((0, 1), (1, 0), 1), ((1, 0), (0, 0), 1), ((0, 0), (1, 1), 1)),
             (((0, 1), (1, 1), 1), ((1, 0), (0, 0), 1)),
             (((0, 1), (1, 1), 1), ((1, 0), (0, 0), 1), ((0, 0), (1, 1), 1))]
        target = set(target)
        self.assertEqual(register, target)
        # Tests second equivalence class
        edges = [((0, 0), (1, 0), 1)]
        graph = create_prime_power_graph(edges, 2, 2)
        class_graph = explore_lc_orbit(graph, verbose=False)
        register = set(tuple(map(tuple, attrs['edges']))
                       for node, attrs in class_graph.node.iteritems())
        target = \
            [(((0, 1), (1, 0), 1),),
             (((0, 1), (1, 0), 1), ((0, 1), (1, 1), 1)),
             (((0, 1), (1, 0), 1), ((0, 1), (1, 1), 1),
              ((1, 0), (0, 0), 1), ((0, 0), (1, 1), 1)),
             (((0, 1), (1, 0), 1), ((1, 0), (0, 0), 1)),
             (((0, 1), (1, 1), 1),),
             (((1, 0), (0, 0), 1),)]
        target = set(target)
        self.assertEqual(register, target)


if __name__ == '__main__':
    unittest.main()
