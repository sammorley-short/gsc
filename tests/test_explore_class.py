# Python packages
import sys
import random
import unittest
import networkx as nx
from abp import GraphState
# Local modules
sys.path.append('..')
from utils import canonical_edge_order
from explore_class import local_complementation


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

    def test_local_complementation(self):
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
            lc_g = local_complementation(g, lc_node)
            lc_gs = to_GraphState(g)
            lc_gs.local_complementation(lc_node)
            # Checks that their edgelists are equal
            lc_g_edges = canonical_edge_order(lc_g.edges())
            lc_gs_edges = canonical_edge_order(lc_gs.edgelist())
            self.assertEqual(lc_g_edges, lc_gs_edges)
            # Performs a second local complementation on same node
            lc_lc_g = local_complementation(lc_g, lc_node)
            lc_lc_gs = to_GraphState(lc_g)
            lc_lc_gs.local_complementation(lc_node)
            # Checks edgelists are equal and also equal to the originals
            lc_lc_g_edges = canonical_edge_order(lc_lc_g.edges())
            lc_lc_gs_edges = canonical_edge_order(lc_lc_gs.edgelist())
            self.assertEqual(lc_lc_g_edges, lc_lc_gs_edges)
            self.assertEqual(lc_lc_g_edges, g_edges)
            self.assertEqual(lc_lc_gs_edges, gs_edges)


if __name__ == '__main__':
    unittest.main()
