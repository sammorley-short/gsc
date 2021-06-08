# Python packages
import pytest
import networkx as nx
import pynauty as pyn
# Local modules
from gsc.get_nauty import find_rep_nodes, hash_graph, canonical_relabel, convert_nx_to_pyn
from gsc.explore_lc_orbit import qubit_LC
from gsc.graph_builders import create_prime_power_graph
from gsc.utils import canonical_edge_order, random_relabel, gen_random_connected_graph

TEST_REPS = 50
LARGE_GRAPH_SIZES = [8, 9, 10]
SMALL_GRAPH_SIZES = [4, 5, 6]


@pytest.mark.parametrize('n_v', LARGE_GRAPH_SIZES)
def test_find_rep_nodes(n_v):
    for _ in range(TEST_REPS):
        g = gen_random_connected_graph(n_v)
        g_equivs = find_rep_nodes(g)
        for _, equiv_nodes in g_equivs.items():
            lc_equiv_graphs = [qubit_LC(g, node) for node in equiv_nodes]
            lc_equiv_hashes = list(set(map(hash_graph, lc_equiv_graphs)))
            assert len(lc_equiv_hashes) == 1


@pytest.mark.parametrize('n_v', LARGE_GRAPH_SIZES)
def test_hash_graph(n_v):
    for _ in range(TEST_REPS):
        g = gen_random_connected_graph(n_v)
        relab_g = random_relabel(g)
        assert hash_graph(g) == hash_graph(relab_g)


@pytest.mark.parametrize('n_v', LARGE_GRAPH_SIZES)
def test_random_relabel(n_v):
    for _ in range(TEST_REPS):
        g = gen_random_connected_graph(n_v)
        relab_g = random_relabel(g)
        assert nx.is_isomorphic(g, relab_g)


@pytest.mark.parametrize('n_v', SMALL_GRAPH_SIZES)
def test_canonical_relabel_random(n_v):
    for _ in range(TEST_REPS):
        g = gen_random_connected_graph(n_v)
        relab_g = random_relabel(g)
        canon_g = canonical_relabel(g)
        canon_relab_g = canonical_relabel(relab_g)
        canon_edges = canonical_edge_order(canon_g.edges())
        canon_relab_edges = canonical_edge_order(canon_relab_g.edges())
        assert canon_edges == canon_relab_edges


@pytest.mark.parametrize('edges_1, edges_2, are_equivalent', [
    # Tests that two C-shaped prime power graphs are equivalent
    (
        [((0, 1), (1, 1), 1), ((1, 0), (0, 0), 1), ((0, 0), (0, 1), 1)],
        [((0, 1), (1, 1), 1), ((1, 0), (0, 0), 1), ((1, 0), (1, 1), 1)],
        True
    ),
    # Tests upper ququart bar and lower ququart bar inequivalent
    (
        [((1, 0), (0, 0), 1)],
        [((0, 1), (1, 1), 1)],
        False
    ),
    # Tests both ququart zigzags are equivalent
    (
        [((0, 1), (1, 1), 1), ((1, 0), (0, 0), 1), ((0, 0), (1, 1), 1)],
        [((0, 1), (1, 1), 1), ((1, 0), (0, 0), 1), ((0, 1), (1, 0), 1)],
        True
    ),
    # Checks two-bar and rotated two-bar ququarts are inequivalent
    (
        [((0, 1), (1, 1), 1), ((0, 0), (1, 0), 1)],
        [((0, 0), (0, 1), 1), ((1, 0), (1, 1), 1)],
        False
    ),
])
def test_prime_power_hash_examples_2_2(edges_1, edges_2, are_equivalent):
    prime, power = 2, 2
    graph_1 = create_prime_power_graph(edges_1, prime, power)
    graph_2 = create_prime_power_graph(edges_2, prime, power)
    if are_equivalent:
        assert hash_graph(graph_1) == hash_graph(graph_2)
    else:
        assert hash_graph(graph_1) != hash_graph(graph_2)


@pytest.mark.parametrize('nx_g_vertices, nx_g_edges, pyn_adj_matrix, expected_pyn_node_map', [
    ([1], [], {0: []}, {0: 1}),
    ([1, 2], [(1, 2)], {0: [1], 1: [0]}, {0: 1, 1: 2}),
    ([2, 4, 6], [(2, 6), (4, 6)], {0: [2], 1: [2], 2: [0, 1]}, {0: 2, 1: 4, 2: 6}),
    ([4, 6, 2], [(2, 6), (4, 6)], {0: [2], 1: [2], 2: [0, 1]}, {0: 2, 1: 4, 2: 6}),
    ([6, 4, 2], [(2, 6), (4, 6)], {0: [2], 1: [2], 2: [0, 1]}, {0: 2, 1: 4, 2: 6}),
])
def test_convert_nx_to_pyn(nx_g_vertices, nx_g_edges, pyn_adj_matrix, expected_pyn_node_map):
    nx_g = nx.Graph()
    nx_g.add_nodes_from(nx_g_vertices)
    nx_g.add_edges_from(nx_g_edges)
    expected_pyn_g = pyn.Graph(len(nx_g_vertices), adjacency_dict=pyn_adj_matrix)
    pyn_g, pyn_g_node_map = convert_nx_to_pyn(nx_g)
    assert pyn_g_node_map == expected_pyn_node_map
    assert pyn_g.adjacency_dict == expected_pyn_g.adjacency_dict
