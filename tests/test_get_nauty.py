# Python packages
import sys
import random
import unittest
import pynauty as pyn
import networkx as nx
# Local modules
sys.path.append('..')
from get_nauty import find_unique_lcs, convert_nx_to_pyn, hash_graph, \
    canonical_relabel
from explore_class import local_complementation
from utils import canonical_edge_order


def gen_random_connected_graph(n, p=0.1):
    """ Generates a Erdos-Renyi G_n,p random graph """
    g = nx.fast_gnp_random_graph(n, p)
    while not nx.is_connected(g):
        g = nx.fast_gnp_random_graph(n, p)
    return g


def random_relabel(graph):
    """ Returns a randomly relabelled graph """
    nodes = list(graph.nodes())
    relab_nodes = list(graph.nodes())
    random.shuffle(relab_nodes)
    relabel = {i_node: o_node for i_node, o_node in zip(nodes, relab_nodes)}
    relabel_graph = nx.relabel_nodes(graph, relabel)
    return relabel_graph


class TestGetNauty(unittest.TestCase):

    def setUp(self):
        pass

    def test_find_unique_lcs(self):
        for _ in range(100):
            g = gen_random_connected_graph(10)
            g_equivs = find_unique_lcs(g)
            for rep_node, equiv_nodes in g_equivs.iteritems():
                lc_equiv_graphs = [local_complementation(g, node)
                                   for node in equiv_nodes]
                lc_equiv_hashes = list(set(map(hash_graph, lc_equiv_graphs)))
                self.assertEqual(len(lc_equiv_hashes), 1)

    def test_hash_graph(self):
        for _ in range(100):
            g = gen_random_connected_graph(10)
            relab_g = random_relabel(g)
            self.assertEqual(hash_graph(g), hash_graph(relab_g))

    def test_random_relabel(self):
        g = gen_random_connected_graph(10)
        relab_g = random_relabel(g)
        self.assertTrue(nx.is_isomorphic(g, relab_g))

    def test_canonical_relabel(self):
        for _ in range(100):
            g = gen_random_connected_graph(4)
            relab_g = random_relabel(g)
            canon_g = canonical_relabel(g)
            canon_relab_g = canonical_relabel(relab_g)
            canon_edges = canonical_edge_order(canon_g.edges())
            canon_relab_edges = canonical_edge_order(canon_relab_g.edges())
            self.assertEqual(canon_edges, canon_relab_edges)


if __name__ == '__main__':
    unittest.main()

    # for _ in range(10):
    #     g = gen_random_connected_graph(4)
    #     relab_g = random_relabel(g)
    #     print g.edges()
    #     print relab_g.edges()
    #     print g.edges() == relab_g.edges()
