# Python packages
from collections import defaultdict
import networkx as nx
import pynauty as pyn
# Local modules
from .utils import (
    compile_maps, invert_dict,
    graph_prime, graph_power, graph_dimension
)
from .qudit_graphs import map_to_qudit_graph


def convert_nx_to_pyn(nx_g, coloring=None):
    """
    Takes a NetworkX graph and outputs a PyNauty graph and a relabelling

    Args:
        nx_g (networkx.Graph): A NetworkX graph.
        coloring (optional, list of lists): The grouping of nodes by color. Should be a list of lists of node
            names where nodes are grouped by colour. For colorless graphs, set to ``None`` (default).

    Returns:
        pynauty.Graph: The equivalent PyNauty graph. Note that this graph will have integer labels.
        dict: The node map from integer nodes in the PyNauty graph to nodes in the NetworkX graph.
    """
    # Check graph is 2D
    if graph_dimension(nx_g) != 2:
        raise RuntimeError('Input graphs must have dimension equal to 2')
    coloring = coloring or []

    # Relabels nodes with integers for compatibility with Pynauty
    nodes = sorted(nx_g.nodes())
    neighs = [sorted(nx_g.adj[node]) for node in nodes]
    to_int_node_map = {n: i for i, n in enumerate(nodes)}
    relabel = to_int_node_map.get
    int_nodes = list(map(relabel, nodes))
    int_neighs = [list(map(relabel, node_neighs)) for node_neighs in neighs]
    int_coloring = [set(map(relabel, colour)) for colour in coloring]

    # Creates Pynauty graph
    graph_adj = dict(zip(int_nodes, int_neighs))
    n_v = len(int_nodes)
    pyn_g = pyn.Graph(n_v, adjacency_dict=graph_adj, vertex_coloring=int_coloring)

    # Finds inverse node labelling
    from_int_node_map = invert_dict(to_int_node_map)

    return pyn_g, from_int_node_map


def hash_graph(graph):
    """ Returns a hash for the graph based on PyNauty's certificate fn """
    if graph_power(graph) > 1:
        return _hash_graph_prime_power_non_trivial(graph)

    if graph_prime(graph) > 2:
        return _hash_graph_prime_power_trivial(graph)

    return _hash_graph_2D(graph)


def _hash_graph_prime_power_non_trivial(graph):
    nx_g_mem, nx_g_mem_coloring = map_to_qudit_graph(graph, partition='member')
    nx_g_fam, nx_g_fam_coloring = map_to_qudit_graph(graph, partition='family')
    pyn_g_mem, _ = convert_nx_to_pyn(nx_g_mem, coloring=nx_g_mem_coloring)
    pyn_g_fam, _ = convert_nx_to_pyn(nx_g_fam, coloring=nx_g_fam_coloring)
    g_mem_hash = hash(pyn.certificate(pyn_g_mem))
    g_fam_hash = hash(pyn.certificate(pyn_g_fam))
    return hash((g_mem_hash, g_fam_hash))


def _hash_graph_prime_power_trivial(graph):
    qudit_g, coloring = map_to_qudit_graph(graph)
    pyn_g, _ = convert_nx_to_pyn(qudit_g, coloring=coloring)
    return hash(pyn.certificate(pyn_g))


def _hash_graph_2D(graph):
    pyn_g, _ = convert_nx_to_pyn(graph)
    return hash(pyn.certificate(pyn_g))


def canonical_relabel(nx_g):
    """ Returns isomorphic graph with canonical relabelling """
    # For edgeless graphs, just return a simple copy
    if not nx_g.edges():
        return nx_g.copy()

    # Generate PyNauty graph and canonically relabel
    pyn_g, pyn_to_nx_relab = convert_nx_to_pyn(nx_g)
    nx_to_pyn_relab = invert_dict(pyn_to_nx_relab)

    # Find out relabelling in integer picture and compile overall relabelling
    pyn_canon_relab = invert_dict(dict(enumerate(pyn.canon_label(pyn_g))))
    nx_canon_relab = compile_maps([nx_to_pyn_relab, pyn_canon_relab, pyn_to_nx_relab])

    return nx.relabel_nodes(nx_g, nx_canon_relab)


def find_rep_nodes(nx_g):
    """
    Takes a NetworkX graph and finds groups of nodes that are equivalent
    up to automorphism
    """
    # Creates PyNauty graph and passes it to PyNauty to get orbits
    partition = 'member' if graph_power(nx_g) > 1 else None
    if graph_dimension(nx_g) != 2:
        nx_g, coloring = map_to_qudit_graph(nx_g, partition)
    else:
        coloring = []
    pyn_g, node_map = convert_nx_to_pyn(nx_g, coloring=coloring)
    _, _, _, orbits, _ = pyn.autgrp(pyn_g)
    # Finds node equivalency dictionary
    node_equivs = defaultdict(list)
    for node, equiv in enumerate(orbits):
        node_equivs[node_map[equiv]].append(node_map[node])
    # Removes any LC's that act trivially on the graph (i.e. act on d=1 nodes)
    node_equivs = {
        node: equivs for node, equivs in node_equivs.items()
        if nx_g.degree[node] > 1
    }
    # If multigraph, returns orbits of nodes in first layer
    if graph_dimension(nx_g) > 2:
        node_equivs = {u: [v for l_v, v in equivs if l_v == 0]
                       for (l_u, u), equivs in node_equivs.items() if l_u == 0}
    return node_equivs
