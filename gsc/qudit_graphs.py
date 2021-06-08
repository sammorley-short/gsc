# Python packages
from math import log
import networkx as nx
# Local modules
from gsc.utils import (
    int_to_bits, copy_graph,
    graph_prime, graph_power
)


def map_to_qudit_graph(nx_wg, partition=None):
    """
    Maps an edge-weighted NetworkX graph to a node-colored NetworkX graph.
    For prime-power graph states, can partition (i.e. color) by member or family.

    Args:
        nx_wg (networkx.Graph): An edge-weighted NetworkX graph.
        partition (str, optional): If the graph is prime-power, then use this argument to specify whether to
            color nodes by 'family' or 'member'. If 'family', nodes within the same qudit are colored the same.
            If 'member' then colors are the same across similarly indexed family members.

    Returns:
        networkx.Graph: A node-colored NetworkX graph.
        list of node lists: The node coloring.
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
    if graph_power(nx_wg) == 1:
        return _get_prime_power_coloring_trivial_power(nx_wg)
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
    if partition == 'family':
        return _get_prime_power_coloring_non_trivial_power_family(nx_wg)

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
