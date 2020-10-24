# Python packages
import networkx as nx
import itertools as it
from random import randint
# Local modules
from gsc.utils import (
    flatten,
    is_prime,
)


def random_connected_graph(n):
    """ Generates a Erdos-Renyi G_{n,m} random graph """
    g = nx.Graph([(0, 1), (2, 3)])
    while not nx.is_connected(g):
        edges = randint(n - 1, n * (n - 1) // 2)
        g = nx.dense_gnm_random_graph(n, edges)
    return g


def linear_graph(l):
    """ Creates a linear graph with 2D coordinates """
    g = nx.Graph()
    edges = [((i, 0), (i + 1, 0)) for i in range(l - 1)]
    g.add_edges_from(edges)
    return g


def square_lattice(n, m, boundary=True):
    """ Creates a square lattice with 2D coordinates """
    mod_n = n + 1 if boundary else n
    mod_m = m + 1 if boundary else m
    g = nx.Graph()
    nodes = it.product(range(n), range(m))
    edges = flatten([[((i, j), ((i + 1) % mod_n, j)),
                      ((i, j), (i, (j + 1) % mod_m))] for i, j in nodes])
    edges = [((u_x, u_y), (v_x, v_y)) for (u_x, u_y), (v_x, v_y) in edges
             if max([u_x, v_x]) < n and max([u_y, v_y]) < m]
    g.add_edges_from(edges)
    return g


def make_crazy(graph, n):
    """ Converts a graph into it's crazily encoded version """
    # Defines groupings of crazy nodes and edges between them
    crazy_nodes = {node: [(node, i) for i in range(n)]
                   for node in graph.nodes()}
    crazy_edges = flatten([it.product(crazy_nodes[u], crazy_nodes[v])
                           for u, v in graph.edges()])
    # Creates graph and adds edges and nodes
    crazy_graph = nx.Graph()
    crazy_graph.add_nodes_from(flatten(crazy_nodes.values()))
    crazy_graph.add_edges_from(crazy_edges)
    # Assigns encoded attribute for to_GraphState
    crazy_graph.encoded = True
    return crazy_graph


def make_ghz_like(graph, n):
    """ Converts a graph into it's GHZ-encoded version """
    # Defines groupings of GHZ-like nodes and edges between them
    ghz_nodes = {node: [(node, i) for i in range(n)]
                 for node in graph.nodes()}
    ghz_edges = [((u, 0), (v, 0)) for u, v in graph.edges()] + \
        [((node, 0), (node, i)) for node in graph.nodes() for i in range(1, n)]
    # Creates graph and adds edges and nodes
    ghz_graph = nx.Graph()
    ghz_graph.add_nodes_from(flatten(ghz_nodes.values()))
    ghz_graph.add_edges_from(ghz_edges)
    # Assigns encoded attribute for to_GraphState
    ghz_graph.encoded = True
    return ghz_graph


def create_prime_graph(w_edges, prime):
    """ Creates a weighted graph representing a prime qudit graph state """
    if not is_prime(prime):
        raise Exception("Graph state must be prime-dimensional")
    us, vs, ws = zip(*w_edges)
    if max(ws) >= prime or max(ws) < 0:
        raise Exception("Weights must be 0 <= w < p ")
    nx_wg = nx.Graph()
    nx_wg.add_weighted_edges_from(w_edges)
    nx_wg.prime = prime
    nx_wg.power = 1
    nx_wg.dimension = prime
    return nx_wg


def create_prime_power_graph(w_edges, prime, power):
    """ Creates a weighted graph representing a prime qudit graph state """
    nx_wg = create_prime_graph(w_edges, prime)
    nx_wg.power, nx_wg.dimension = power, prime ** power
    fam_labels = list(set([n for n, i in nx_wg.nodes()]))
    nx_wg.families = len(fam_labels)
    fam_nodes = [(n, i) for n in fam_labels for i in range(power)]
    # Adds any nodes that weren't in the edge list
    nx_wg.add_nodes_from(fam_nodes)
    return nx_wg


def from_MDS_code(A, prime, power):
    """ Creates a graph-state representation of an MDS AME state """
    A = np.array(A)
    A_t = np.transpose(A)
    k, nmk = A.shape
    t_l_zero = np.zeros((k, k), dtype=int)
    b_r_zero = np.zeros((nmk, nmk), dtype=int)
    top = np.concatenate((t_l_zero, A), axis=1)
    bot = np.concatenate((A_t, b_r_zero), axis=1)
    adj_mat = np.concatenate((top, bot), axis=0)
    graph = nx.from_numpy_array(adj_mat)
    graph.prime = prime
    graph.power = power
    graph.dimension = prime ** power
    return graph
