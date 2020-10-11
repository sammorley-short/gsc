# Python packages
import random
import networkx as nx
# Local modules
from gsc.get_nauty import find_rep_nodes, hash_graph, canonical_relabel
from gsc.explore_lc_orbit import qubit_LC
from gsc.graph_builders import create_prime_power_graph
from gsc.utils import canonical_edge_order


def test_find_rep_nodes():
    for _ in range(100):
        g = gen_random_connected_graph(10)
        g_equivs = find_rep_nodes(g)
        for rep_node, equiv_nodes in g_equivs.iteritems():
            lc_equiv_graphs = [qubit_LC(g, node)
                               for node in equiv_nodes]
            lc_equiv_hashes = list(set(map(hash_graph, lc_equiv_graphs)))
            assert len(lc_equiv_hashes) == 1


def test_hash_graph():
    for _ in range(100):
        g = gen_random_connected_graph(10)
        relab_g = random_relabel(g)
        assert hash_graph(g) == hash_graph(relab_g)


def test_random_relabel():
    g = gen_random_connected_graph(10)
    relab_g = random_relabel(g)
    assert nx.is_isomorphic(g, relab_g)


def test_canonical_relabel():
    for _ in range(100):
        g = gen_random_connected_graph(4)
        relab_g = random_relabel(g)
        canon_g = canonical_relabel(g)
        canon_relab_g = canonical_relabel(relab_g)
        canon_edges = canonical_edge_order(canon_g.edges())
        canon_relab_edges = canonical_edge_order(canon_relab_g.edges())
        assert canon_edges == canon_relab_edges


def test_prime_power_hash_examples():
    # Tests that two C-shaped prime power graphs are equivalent
    prime, power = 2, 2
    e1 = [((0, 1), (1, 1), 1), ((1, 0), (0, 0), 1), ((0, 0), (0, 1), 1)]
    e2 = [((0, 1), (1, 1), 1), ((1, 0), (0, 0), 1), ((1, 0), (1, 1), 1)]
    g1 = create_prime_power_graph(e1, prime, power)
    g2 = create_prime_power_graph(e2, prime, power)
    assert hash_graph(g1) == hash_graph(g2)

    # Tests upper ququart bar and lower ququart bar inequivalent
    e1 = [((1, 0), (0, 0), 1)]
    e2 = [((0, 1), (1, 1), 1)]
    g1 = create_prime_power_graph(e1, prime, power)
    g2 = create_prime_power_graph(e2, prime, power)
    assert hash_graph(g1) == hash_graph(g2)

    # Tests both ququart zigzags are equivalent
    e1 = [((0, 1), (1, 1), 1), ((1, 0), (0, 0), 1), ((0, 0), (1, 1), 1)]
    e2 = [((0, 1), (1, 1), 1), ((1, 0), (0, 0), 1), ((0, 1), (1, 0), 1)]
    g1 = create_prime_power_graph(e1, prime, power)
    g2 = create_prime_power_graph(e2, prime, power)
    assert hash_graph(g1) == hash_graph(g2)

    # Checks two-bar and rotated two-bar ququarts are inequivalent
    e1 = [((0, 1), (1, 1), 1), ((0, 0), (1, 0), 1)]
    e2 = [((0, 0), (0, 1), 1), ((1, 0), (1, 1), 1)]
    g1 = create_prime_power_graph(e1, prime, power)
    g2 = create_prime_power_graph(e2, prime, power)
    assert hash_graph(g1) == hash_graph(g2)


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
