# Python packages
from math import log
import networkx as nx
import pynauty as pyn
from collections import defaultdict
# Local modules
from gsc.utils import (
    int_to_bits, copy_graph, compile_maps, invert_dict,
    graph_prime, graph_power, graph_dimension
)


def qudit_graph_map(nx_wg, partition=None):
    """
    Maps an edge-weighted NetworkX graph to a node-colored NetworkX graph.
    For prime-power graph states, can partition (i.e. color) by member or family.
    """

    # Family-partitioned graphs need to be extended to allow for exchanges of entire families
    if partition == 'family' and graph_power(nx_wg) > 1:
        # If we're modifying, copy graph to avoid in-place edits of original
        nx_wg = copy_graph(nx_wg)
        _extend_prime_power_graph_non_trivial_power(nx_wg)

    # Generate qudit graph and it's coloring
    nx_cg = _create_qudit_graph(nx_wg)
    coloring = _get_prime_power_coloring(nx_wg, partition)

    # Creates layered graph with vertical edges
    return nx_cg, coloring


def _extend_prime_power_graph_non_trivial_power(nx_wg):
    """
    Adds extra nodes to represent exchangeable colours
    (see page 60 of nauty user guide v26)
    """
    for family in range(nx_wg.families):
        node = (family, nx_wg.power)
        nx_wg.add_node(node)
        nx_wg.add_weighted_edges_from([
            (node, (family, i), 1)
            for i in range(nx_wg.power)
        ])


def _create_qudit_graph(nx_wg):
    layers = _get_edge_weight_layers(nx_wg)
    nx_cg = nx.Graph()
    nx_cg.add_nodes_from([
        (layer, i)
        for layer in layers
        for i in nx_wg.nodes()
    ])
    nx_cg.add_edges_from([
        ((layer, i), (layer + 1, i))
        for i in nx_wg.nodes()
        for layer in layers[:-1]
    ])

    # Add edges within layers
    for vertex_i, vertex_j, weight in nx_wg.edges.data('weight'):
        # Gets binary rep. of weight, padded with zeros (written L to R)
        bin_w = int_to_bits(weight)[::-1]
        bin_w = bin_w + (len(layers) - len(bin_w)) * [0]

        # Converts binary weight to list of layers and adds edges to graph
        edge_layers = [layer for layer, i in enumerate(bin_w) if i]
        edges = [((layer, vertex_i), (layer, vertex_j)) for layer in edge_layers]
        nx_cg.add_edges_from(edges)

    return nx_cg


def _get_prime_power_coloring(nx_wg, partition):
    # If node is prime power, applies colouring across same member-nodes
    if graph_prime(nx_wg) == 1:
        return _get_prime_power_coloring_trivial_power(nx_wg)
    else:
        return _get_prime_power_coloring_non_trivial_power(nx_wg, partition)


def _get_prime_power_coloring_trivial_power(nx_wg):
    return [
        [
            (layer, n)
            for n in nx_wg.nodes()
        ]
        for layer in _get_edge_weight_layers(nx_wg)
    ]


def _get_edge_weight_layers(nx_wg):
    """
    Gets the list of layer indices needed to represent edge weights/colours.
    See page 60 of nauty user guide v27 for more details.
    """
    _, _, weights = zip(*nx_wg.edges.data('weight'))
    n_layers = int(log(max(weights), 2)) + 1
    return range(n_layers)


def _get_prime_power_coloring_non_trivial_power(nx_wg, partition):
    # Partition based on which member of the family node is
    if partition == 'member':
        return _get_prime_power_coloring_non_trivial_power_member(nx_wg)
    # Partitions based on which family => must make colourings equiv.
    elif partition == 'family':
        return _get_prime_power_coloring_non_trivial_power_family(nx_wg)
    else:
        raise Exception(f"Unknown partition scheme {partition} provided; choose from 'member' or 'family'.")


def _get_prime_power_coloring_non_trivial_power_member(nx_wg):
    return [
        [
            (layer, (family, i))
            for family in range(nx_wg.families)
        ]
        for layer in _get_edge_weight_layers(nx_wg)
        for i in range(nx_wg.power)
    ]


def _get_prime_power_coloring_non_trivial_power_family(nx_wg):
    # N.B. recall that m > 1 family-coloured graphs will have been extended with new nodes to allow for
    # exchangeable families.
    layers = _get_edge_weight_layers(nx_wg)
    power, families = nx_wg.power, nx_wg.families

    old_node_colorings = [
        [
            (layer, (family, i))
            for family in range(families)
            for i in range(power)
        ]
        for layer in layers
    ]
    new_node_colorings = [
        [
            (layer, (family, power))
            for family in range(families)
        ]
        for layer in layers
    ]
    return old_node_colorings + new_node_colorings


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
        nx_g_mem, nx_g_mem_coloring = qudit_graph_map(graph, partition='member')
        nx_g_fam, nx_g_fam_coloring = qudit_graph_map(graph, partition='family')
        pyn_g_mem, _ = convert_nx_to_pyn(nx_g_mem, coloring=nx_g_mem_coloring)
        pyn_g_fam, _ = convert_nx_to_pyn(nx_g_fam, coloring=nx_g_fam_coloring)
        g_mem_hash = hash(pyn.certificate(pyn_g_mem))
        g_fam_hash = hash(pyn.certificate(pyn_g_fam))
        g_hash = hash((g_mem_hash, g_fam_hash))
    else:
        pyn_g, _ = convert_nx_to_pyn(graph)
        g_hash = hash(pyn.certificate(pyn_g))
    return g_hash


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
    partition = 'member' if nx_g.__dict__.get('power', 1) > 1 else None
    if nx_g.__dict__.get('dimension', 2) != 2:
        nx_g, coloring = qudit_graph_map(nx_g, partition)
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
    if nx_g.__dict__.get('dimension', 2) > 2:
        node_equivs = {u: [v for l_v, v in equivs if l_v == 0]
                       for (l_u, u), equivs in node_equivs.items() if l_u == 0}
    return node_equivs
