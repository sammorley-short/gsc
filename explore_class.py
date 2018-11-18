# Python packages
import os
import json
import pynauty as pyn
import networkx as nx
import itertools as it
from copy import deepcopy
from networkx.readwrite import json_graph
from pprint import pprint
# Local modules
from utils import canonical_edge_order
from get_nauty import find_unique_lcs, hash_graph


def local_complementation(graph, node, copy=True):
    """ Returns the graph for local complementation applied to node """
    neighs = graph.neighbors(node)
    neigh_k_edges = it.combinations(neighs, 2)
    lc_graph = deepcopy(graph) if copy else graph
    for u, v in neigh_k_edges:
        if lc_graph.has_edge(u, v):
            lc_graph.remove_edge(u, v)
        else:
            lc_graph.add_edge(u, v)
    return lc_graph


def edge_lc(graph, edge):
    u, v = edge
    edge_lc_graph = apply_lcs(graph, [u, v, u])
    return edge_lc_graph


def apply_lcs(graph, nodes):
    lc_graph = deepcopy(graph)
    for node in nodes:
        lc_graph = local_complementation(lc_graph, node, copy=False)
    return lc_graph


def init_EC_database_dir(directory='EC_database'):
    """ Initialises the equivalence class database directory """
    if not os.path.exists(directory):
        os.makedir(directory)


def recursive_LC_search(class_graph, graph_label, node_equivs):
    """ Recursive search for exploring the isomorphic LC space of graphs """
    graph = class_graph.node[graph_label]['nx_graph']
    for rep_node, equiv_nodes in node_equivs.iteritems():
        lc_graph = local_complementation(graph, rep_node)
        lc_edges = canonical_edge_order(lc_graph.edges())
        try:
            lc_label = class_graph.member_hash_table[lc_edges]
            class_graph.add_edge(graph_label, lc_label,
                                 equivs=equiv_nodes)
            continue
        except KeyError:
            lc_label = max(class_graph.nodes()) + 1
            lc_graph_hash = hash_graph(lc_graph)
            class_graph.add_node(lc_label, nx_graph=lc_graph,
                                 edges=lc_edges, hash=lc_graph_hash)
            class_graph.add_edge(graph_label, lc_label,
                                 equivs=equiv_nodes)
            class_graph.member_hash_table.update({lc_edges: lc_label})
            lc_equivs = find_unique_lcs(lc_graph)
            recursive_LC_search(class_graph, lc_label, lc_equivs)
    return class_graph


def int_relabel_graph(graph):
    """
        Relabels graphs with tuple node names to int node names.
        Returns relabelled graph and node mapping applied.
    """
    int_labels = {node: i for i, node in enumerate(graph.nodes())}
    int_graph = nx.relabel_nodes(graph, int_labels)
    return int_graph, int_labels


def explore_LC_isomorphic_orbit(init_graph):
    """ Explores the LC equivalence class orbit up to isomorphism """
    if not nx.is_connected(init_graph):
        raise TypeError("Initial graph must be connected.")
    init_graph_equivs = find_unique_lcs(init_graph)
    init_edges = canonical_edge_order(init_graph.edges())
    init_hash = hash_graph(init_graph)
    class_graph = nx.Graph()
    class_graph.add_node(0, nx_graph=init_graph, edges=init_edges,
                         hash=init_hash)
    class_graph.member_hash_table = {init_edges: 0}
    class_graph = recursive_LC_search(class_graph, 0, init_graph_equivs)
    return class_graph


def export_class_graph(class_graph, filename):
    """ Exports class graph to JSON file """
    for node, attrs in class_graph.node.iteritems():
        class_graph.node[node].pop('nx_graph')
    class_graph_data = json_graph.adjacency_data(class_graph)
    with open(filename, 'w') as fp:
        json.dump(class_graph_data, fp)
    return class_graph_data


if __name__ == '__main__':
    n = 4
    init_edges = [(i, (i + 1) % n) for i in range(n)]
    init_graph = nx.Graph(init_edges)
    print init_graph.edges()
    print
    class_graph = explore_LC_isomorphic_orbit(init_graph)
    pprint(dict(class_graph.node))
