# Python modules
import networkx as nx
import itertools as it
# Local modules
from gsc.utils import powerset


def gen_psuedo_graph_edge_map(prime, power):
    """
    Creates mapping from p^m psuedo graph edges to prime-dimensional graph
    state edges (i.e. set of all weighted balanced bipartite graphs between
    two m-families).
    """
    edges = ([(u, v, w) for w in range(1, prime)] + [()]
             for u, v in it.product(range(power), repeat=2))
    edge_sets = ([edge for edge in edge_set if edge != ()]
                 for edge_set in it.product(*edges))
    c_map = {i: edges for i, edges in enumerate(sorted(edge_sets, key=len))}
    return c_map


def create_psuedo_graph(c_edges, prime, power, c_map):
    """
    Initialises prime power psuedo graph state.
    c_edges is list of coloured/weighted edges that are associated to bipartite
    prime qudit graphs via c_map
    """
    # Initialises directed weighted psuedo graph and adds attributes
    nx_wg = nx.DiGraph()
    # Adds psuedo edges and real prime edges as edge attribute
    nx_wg.add_weighted_edges_from(c_edges)
    nx_wg.remove_edges_from([(u, v) for u, v, c in c_edges if c == 0])
    # Add psuedo graph
    nx_wg.prime, nx_wg.power, nx_wg.dimension = prime, power, prime ** power
    nx_wg.c_map = c_map
    return nx_wg


def psuedo_to_real(psu_g):
    """ Turns psuedo-graph into real prime version """
    prime = psu_g.prime
    power = psu_g.power
    c_map = psu_g.c_map
    real_g = nx.Graph()
    nodes = [(n, i) for n in psu_g.nodes() for i in range(power)]
    real_g.add_nodes_from(nodes)
    edges = [((u, i), (v, j), w) for u, v, c in psu_g.edges(data='weight')
             for i, j, w in c_map[c]]
    real_g.add_weighted_edges_from(edges)
    real_g.families = len(psu_g.nodes())
    real_g.prime, real_g.power, real_g.dimension = prime, power, prime ** power
    return real_g


def real_graph_to_psu_edges(real_g, c_map, psu_edge_index):
    """
    Converts real graph to list of psuedo edges ordered by the pseudo edge
    index
    """
    # Creates inverse colour map for psuedo edges
    inv_c_map = {tuple(sorted(value)): key for key, value in c_map.items()}
    # Gets configuration of real edges per psuedo edge
    psu_edges = {edge: [] for edge in psu_edge_index}
    m = real_g.power
    for u, v in psu_edge_index:
        # Gets edges of inter-family bipartite graph
        bpg_edges = [(i, j, real_g[(u, i)][(v, j)]['weight'])
                     for i, j in it.product(range(m), repeat=2)
                     if real_g.has_edge((u, i), (v, j))]
        if bpg_edges != []:
            psu_edges[(u, v)] = bpg_edges
    # Creates weighted directed graph from psuedo edges
    psu_edges = [(u, v, inv_c_map[tuple(sorted(psu_edges[(u, v)]))])
                 for u, v in psu_edge_index]
    return psu_edges


def real_to_psuedo(real_g, c_map, psu_edge_index=None):
    """ Converts a real graph back to a psuedo graph given a mapping func """
    # Finds value for every possible psuedo edge
    psu_nodes = list(set([u for u, i in real_g.nodes()]))
    all_psu_edges = psu_edge_index if psu_edge_index \
        else it.combinations(psu_nodes, 2)
    psu_edges = real_graph_to_psu_edges(real_g, c_map, all_psu_edges)
    psu_g = nx.DiGraph()
    psu_g.prime, psu_g.power = real_g.prime, psu_g.power
    psu_g.dimension, psu_g.c_map = real_g.dimension, c_map
    psu_g.add_weighted_edges_from(psu_edges)
    return psu_g


if __name__ == '__main__':
    prime, power = 3, 1
    print gen_psuedo_graph_edge_map(prime, power)
