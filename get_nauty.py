# Python packages
import networkx as nx
import pynauty as pyn
from pprint import pprint
from collections import defaultdict
# Local modules
from utils import flatten


def convert_nx_to_pyn(nx_g):
    """ Takes a NetworkX graph and outputs a PyNauty graph """
    nodes, neighs = zip(*nx_g.adjacency())
    graph_adj = {node: neighs.keys() for node, neighs in zip(nodes, neighs)}
    # print graph_adj
    n_v = len(graph_adj)
    pyn_g = pyn.Graph(n_v, directed=False, adjacency_dict=graph_adj)
    return pyn_g


def hash_graph(graph):
    """ Returns a hash for the graph based on PyNauty's certificate fn """
    pyn_g = convert_nx_to_pyn(graph)
    g_cert = pyn.certificate(pyn_g)
    return hash(g_cert)


def canonical_relabel(nx_g):
    """ Returns isomorphic graph with canonical relabelling """
    nodes, neighs = zip(*nx_g.adjacency())
    pyn_g = convert_nx_to_pyn(nx_g)
    canon_lab = pyn.canon_label(pyn_g)
    canon_relab = {o_node: i_node for i_node, o_node in zip(nodes, canon_lab)}
    nx_g_canon = nx.relabel_nodes(nx_g, canon_relab)
    return nx_g_canon


def find_automorphisms(nx_g):
    """
        Takes a NetworkX graph and returns the automorphism generators
    """
    nodes, neighs = zip(*nx_g.adjacency())
    pyn_g = convert_nx_to_pyn(nx_g)
    gens, grp1_size, grp2_size, orbits, n_orbits = pyn.autgrp(pyn_g)
    am_node_maps = [{i_node: o_node for i_node, o_node in zip(nodes, gen)}
                    for gen in gens]
    return nx_g, am_node_maps


def find_unique_lcs(g):
    """
        Takes a NetworkX graph and finds groups of nodes that are equivalent
        up to automorphism
    """
    g, am_node_maps = find_automorphisms(g)
    aut_edges = flatten([am_node_map.items() for am_node_map in am_node_maps])
    aut_g = nx.Graph(aut_edges)
    equiv_groups = [equiv for equiv in map(list, list(
        nx.connected_components(aut_g))) if len(equiv) > 1]
    all_nodes = g.nodes()
    node_equivs = {}
    for equiv_group in equiv_groups:
        rep_node = equiv_group[0]
        node_equivs.update({rep_node: equiv_group})
        all_nodes = [node for node in all_nodes if node not in equiv_group]
    node_equivs.update({node: [node] for node in all_nodes})
    # Removes any LC's that act trivially on the graph (i.e. act on d=1 nodes)
    node_equivs = {node: equivs for node, equivs in node_equivs.iteritems()
                   if g.degree(node) > 1}
    return node_equivs


# TODO: Comment functions
# TODO: Package up hacked nauty for distribution


if __name__ == '__main__':
    n = 12
    edge_sets = [[(i, (i + 1) % n) for i in range(n)],
                 [(i, (i + 1) % n) for i in range(n)] + [(1, 5)],
                 [(i, (i + 1)) for i in range(n)] + [(2, 7)]]
    for edges in edge_sets:
        g = nx.Graph(edges)
        print "Edges:"
        pprint(list(g.edges()))
        # node_equivs = find_unique_lcs(g)
        # print "Node LC equivalences:", node_equivs
        # print hash_graph(g)
        relab_g = canonical_relabel(g)
        print "Relabelled edges:"
        pprint(list(relab_g.edges()))
        rr_g = canonical_relabel(relab_g)
        print "Relabelled edges:"
        pprint(list(rr_g.edges()))
        print

    # g = nx.Graph([(0, 2), (0, 3), (0, 5), (1, 4), (1, 5), (3, 4), (3, 5)])
    # g, g_ams = find_automorphisms(g)
    # print g_ams
