# Python packages
from math import log
import networkx as nx
import pynauty as pyn
from collections import defaultdict
# Local modules
from gsc.utils import int_to_bits, copy_graph, compile_maps, invert_dict


def qudit_graph_map(nx_wg, partition=None):
    """
    Maps an edge-weighted NX graph to a node-colored NX graph.
    For prime-power graph states, can colour by member or family.
    """
    # Gets list of all nodes by layer
    us, vs, weights = zip(*nx_wg.edges.data('weight'))
    n_layers = int(log(max(weights), 2)) + 1
    layers = range(n_layers)

    # If node is prime power, applies colouring across same member-nodes
    if nx_wg.__dict__.get('power', 1) > 1:
        _, m, f = nx_wg.prime, nx_wg.power, nx_wg.families

        # Partitions based on which member of the family node is
        if partition == 'member':
            coloring = [
                [
                    (l, (n, i))
                    for n in range(f)
                ]
                for l in layers
                for i in range(m)
            ]

        # Partitions based on which family => must make colourings equiv.
        elif partition == 'family':

            # Adds extra nodes to represent exchangeable colours
            # (see page 60 of nauty user guide v26)
            nx_wg = copy_graph(nx_wg)
            for u in range(f):
                node = (u, m)
                nx_wg.add_node(node)
                nx_wg.add_weighted_edges_from([
                    (node, (u, i), 1)
                    for i in range(m)
                ])
            coloring = [
                [
                    (l, (n, i))
                    for n in range(f)
                    for i in range(m)
                ]
                for l in layers
            ] + [
                [
                    (l, (n, m))
                    for n in range(f)
                ]
                for l in layers
            ]

        else:
            raise Exception("Unknown colour scheme provided")

    else:
        coloring = [
            [(l, n) for n in nx_wg.nodes()]
            for l in layers
        ]

    # Creates layered graph with vertical edges
    nx_cg = nx.Graph()
    v_nodes = [
        (l, n)
        for l in layers
        for n in nx_wg.nodes()
    ]
    nx_cg.add_nodes_from(v_nodes)
    v_edges = [
        ((l, n), (l + 1, n))
        for n in nx_wg.nodes()
        for l in layers[:-1]
    ]
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
    if nx_g.__dict__.get('dimension', 2) != 2:
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
    if graph.__dict__.get('power', 1) > 1:
        pyn_g_mem, _ = convert_nx_to_pyn(graph, partition='member')
        pyn_g_fam, _ = convert_nx_to_pyn(graph, partition='family')
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
