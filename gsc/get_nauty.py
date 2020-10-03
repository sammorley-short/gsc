# Python packages
from math import log
import networkx as nx
import pynauty as pyn
import zlib
from collections import defaultdict
# Local modules
from gsc.utils import int_to_bits, copy_graph


def qudit_graph_map(nx_wg, partition=None):
    """
    Maps an edge-weighted NX graph to a node-colored NX graph.
    For prime-power graph states, can colour by member or family.
    """
    # Gets list of all nodes by layer
    us, vs, weights = list(zip(*nx_wg.edges.data('weight')))
    n_layers = int(log(max(weights), 2)) + 1
    layers = list(range(n_layers))
    # If node is prime power, applies colouring across same member-nodes
    if nx_wg.__dict__.get('power', 1) > 1:
        _, m, f = nx_wg.prime, nx_wg.power, nx_wg.families
        # Partitions based on which member of the family node is
        if partition == 'member':
            coloring = [[(l, (n, i)) for n in range(f)]
                        for l in layers for i in range(m)]
        # Partitions based on which family => must make colourings equiv.
        elif partition == 'family':
            # Adds extra nodes to represent exchangeable colours
            # (see page 60 of nauty user guide v26)
            nx_wg = copy_graph(nx_wg)
            for u in range(f):
                node = (u, m)
                nx_wg.add_node(node)
                nx_wg.add_weighted_edges_from([(node, (u, i), 1)
                                               for i in range(m)])
            coloring = [[(l, (n, i)) for n in range(f) for i in range(m)]
                        for l in layers] + \
                       [[(l, (n, m)) for n in range(f)] for l in layers]
        else:
            raise Exception("Unknown colour scheme provided")
    else:
        coloring = [[(l, n) for n in nx_wg.nodes()] for l in layers]
    # Creates layered graph with vertical edges
    nx_cg = nx.Graph()
    v_nodes = [(l, n) for l in layers for n in nx_wg.nodes()]
    nx_cg.add_nodes_from(v_nodes)
    v_edges = [((l, n), (l + 1, n))
               for n in nx_wg.nodes() for l in layers[:-1]]
    nx_cg.add_edges_from(v_edges)
    # Add edges within layers
    for u, v, w in nx_wg.edges.data('weight'):
        # Gets binary rep. of weight, padded with zeros (written L to R)
        bin_w = int_to_bits(w)[::-1]
        bin_w = bin_w + (n_layers - len(bin_w)) * [0]
        # Converts binary weight to list of layers and adds edges to graph
        edge_layers = [l for l, b in enumerate(bin_w) if b]
        edges = [((l, u), (l, v)) for l in edge_layers]
        nx_cg.add_edges_from(edges)
    return nx_cg, coloring


def convert_nx_to_pyn(nx_g, partition=None):
    """
        Takes a NetworkX graph and outputs a PyNauty graph.
        If graph has dimension > 2, converts into layered coloured graph
    """
    # If graph represents nD qudit graph, map to coloured layer graph
    coloring = []
    if nx_g.__dict__.get('dimension', 2) > 2:
        nx_g, coloring = qudit_graph_map(nx_g, partition)
    # Relabels nodes with integers for compatibility with Pynauty
    nodes, neighs = list(zip(*nx_g.adjacency()))
    to_int_node_map = {n: i for i, n in enumerate(nodes)}
    relabel = to_int_node_map.get
    nodes = list(map(relabel, nodes))
    neighs = [list(map(relabel, list(node_neighs.keys()))) for node_neighs in neighs]
    coloring = [set(map(relabel, colour)) for colour in coloring]
    # Creates Pynauty graph
    graph_adj = {node: node_neighs for node, node_neighs in zip(nodes, neighs)}
    n_v = len(graph_adj)
    pyn_g = pyn.Graph(n_v, directed=False, adjacency_dict=graph_adj,
                      vertex_coloring=coloring)
    # Finds inverse node labelling
    from_int_node_map = {i: n for n, i in list(to_int_node_map.items())}
    return pyn_g, from_int_node_map


def hash_graph(graph):
    """ Returns a hash for the graph based on PyNauty's certificate fn """
    if graph.__dict__.get('power', 1) > 1:
        pyn_g_mem, _ = convert_nx_to_pyn(graph, partition='member')
        pyn_g_fam, _ = convert_nx_to_pyn(graph, partition='family')

        g_mem_hash = zlib.crc32(pyn.certificate(pyn_g_mem))
        g_fam_hash = zlib.crc32(pyn.certificate(pyn_g_fam))
        g_hash = zlib.crc32((g_mem_hash, g_fam_hash))
    else:
        pyn_g, _ = convert_nx_to_pyn(graph)
        g_hash = zlib.crc32(pyn.certificate(pyn_g))
    return g_hash


def canonical_relabel(nx_g):
    """ Returns isomorphic graph with canonical relabelling

    Note: Function is broken and not actually used either.

    """
    nodes, neighs = list(zip(*nx_g.adjacency()))
    pyn_g, node_map = convert_nx_to_pyn(nx_g)
    canon_lab = pyn.canon_label(pyn_g)
    canon_relab = {node_map[o_node]: i_node for i_node, o_node
                   in zip(nodes, canon_lab)}
    nx_g_canon = nx.relabel_nodes(nx_g, canon_relab)
    return nx_g_canon


def find_rep_nodes(nx_g):
    """
    Takes a NetworkX graph and finds groups of nodes that are equivalent
    up to automorphism
    """
    # Creates PyNauty graph and passes it to PyNauty to get orbits
    partition = 'member' if nx_g.__dict__.get('power', 1) > 1 else None
    pyn_g, node_map = convert_nx_to_pyn(nx_g, partition=partition)
    _, _, _, orbits, _ = pyn.autgrp(pyn_g)
    # Finds node equivalency dictionary
    node_equivs = defaultdict(list)
    for node, equiv in enumerate(orbits):
        node_equivs[node_map[equiv]].append(node_map[node])

    # If multigraph, returns orbits of nodes in first layer
    if nx_g.__dict__.get('dimension', 2) > 2:
        node_equivs = {u: [v for l_v, v in equivs if l_v == 0]
                       for (l_u, u), equivs in list(node_equivs.items()) if
                       l_u == 0}
    else:
        node_equivs = {node: equivs for node, equivs in node_equivs.items()
                       if nx_g.degree(node) > 1}

    return node_equivs
