# Python packages
import os
import sys
import csv
import json
import pynauty as pyn
import networkx as nx
import itertools as it
from pprint import pprint
from networkx.readwrite import json_graph
# Local modules
from gsc.utils import *
from gsc.get_nauty import find_rep_nodes, hash_graph


def init_EC_database_dir(directory='EC_database'):
    """ Initialises the equivalence class database directory """
    if not os.path.exists(directory):
        os.makedirs(directory)


def qubit_LC(graph, node, copy=True):
    """ Returns the graph for local complementation applied to node """
    neighs = graph.neighbors(node)
    neigh_k_edges = it.combinations(neighs, 2)
    lc_graph = copy_graph(graph) if copy else graph
    for u, v in neigh_k_edges:
        if lc_graph.has_edge(u, v):
            lc_graph.remove_edge(u, v)
        else:
            lc_graph.add_edge(u, v, weight=1)
    return lc_graph


def apply_qubit_LCs(graph, nodes):
    """ Applies a sequence of local complementations """
    lc_graph = copy_graph(graph)
    for node in nodes:
        lc_graph = qubit_LC(lc_graph, node, copy=False)
    return lc_graph


def edge_LC(graph, edge):
    """ Applies edge-local complementation """
    u, v = edge
    edge_lc_graph = apply_qubit_LCs(graph, [u, v, u])
    return edge_lc_graph


def prime_qudit_LC(graph, node, a, copy=True):
    """
    Returns the graph for generalised local complementation applied to
    node n with weight a
    """
    p = graph.prime
    neighs = list(graph.neighbors(node))
    neigh_k_edges = it.combinations(neighs, 2)
    lc_graph = copy_graph(graph) if copy else graph
    for u, v in neigh_k_edges:
        nu_weight = lc_graph[node][u]['weight']
        nv_weight = lc_graph[node][v]['weight']
        if lc_graph.has_edge(u, v):
            uv_weight = lc_graph[u][v]['weight']
            new_weight = (uv_weight + a * nv_weight * nu_weight) % p
            lc_graph[u][v]['weight'] = new_weight
        else:
            weight = (a * nv_weight * nu_weight) % p
            lc_graph.add_edge(u, v, weight=weight)
        if lc_graph[u][v]['weight'] == 0:
            lc_graph.remove_edge(u, v)
    return lc_graph


def prime_qudit_EM(graph, node, b, copy=True):
    """
    Returns the graph for edge multiplication operation applied on node n
    with weight b
    """
    em_graph = copy_graph(graph) if copy else graph
    p = graph.prime
    neighs = em_graph.neighbors(node)
    for u in neighs:
        if em_graph.has_edge(node, u):
            nv_weight = em_graph[node][u]['weight']
            new_weight = (b * nv_weight) % p
            em_graph[node][u]['weight'] = new_weight
        if em_graph[node][u]['weight'] == 0:
            em_graph.remove_edge(node, u)
    return em_graph


def make_LC_a(a):
    def lc_a(graph, node):
        return prime_qudit_LC(graph, node, a)
    return lc_a


def make_EM_b(b):
    def em_b(graph, node):
        return prime_qudit_EM(graph, node, b)
    return em_b


def prime_power_qudit_CC(graph, node, t, a, copy=True):
    """ Returns graph after controlled complementation """
    # Creates new graph if needed
    new_graph = copy_graph(graph) if copy else graph
    n, c = node
    # Adds edge between control and target node
    if c != t:
        new_graph.add_edge((n, c), (n, t), weight=1)
    # Applies LC to control
    new_graph = prime_qudit_LC(new_graph, (n, c), a, copy=False)
    # Removes any intra-family edges
    family_edges = [((u, i), (v, j)) for (u, i), (v, j) in new_graph.edges()
                    if u == v]
    new_graph.remove_edges_from(family_edges)
    return new_graph


def make_pp_CC_a(a, t):
    """ Controlled complementation function generator """
    def cc_a(graph, node):
        return prime_power_qudit_CC(graph, node, t, a)
    return cc_a


def queued_orbit_search(init_graph, local_ops, save_edges, verbose):
    # Initialises class graph with init_graph
    init_edges = list(init_graph.edges())
    init_hash = hash_graph(init_graph)
    class_graph = nx.Graph()
    class_graph.add_node(0, nx_graph=init_graph,
                         edges=init_edges, hash=init_hash)
    class_graph.member_hash_table = {init_hash: 0}
    # Loops over queue members until empty
    queue = [0]
    visited = 0
    while queue:
        # Prints live count of explored/known
        if verbose:
            out = \
                str(visited) + '/' + str(len(queue) + visited) + ' visited (' \
                + str(int(100 * float(visited)/(len(queue) + visited))) + '%)'
            sys.stdout.write('%s\r' % out)
            sys.stdout.flush()
        visited += 1
        # Gets next graph on queue and finds representative nodes
        graph_label = queue.pop()
        graph = class_graph.node[graph_label]['nx_graph']
        node_equivs = find_rep_nodes(graph)
        # Applies set of local ops to each representative node
        for rep_node, equiv_nodes in node_equivs.iteritems():
            for op_label, local_op in local_ops:
                new_graph = local_op(graph, rep_node)
                new_edges = list(new_graph.edges())
                new_hash = hash_graph(new_graph)
                # Checks new graph is difference to original
                if sorted(new_graph.edges(data='weight')) == \
                        sorted(graph.edges(data='weight')):
                    continue
                # If different, tries to find new graph in class
                try:
                    old_label = class_graph.member_hash_table[new_hash]
                    if save_edges:
                        # If new edge adds new edge between members
                        if not class_graph.has_edge(graph_label, old_label):
                            class_graph.add_edge(graph_label, old_label,
                                                 equivs=[equiv_nodes],
                                                 ops=[op_label])
                        # Else adds any new local ops to edge label
                        elif op_label not in \
                                class_graph[graph_label][old_label]['ops']:
                            class_graph[graph_label][old_label]['ops']\
                                .append(op_label)
                            class_graph[graph_label][old_label]['equivs']\
                                .append(equiv_nodes)
                    continue
                # If not in class, creates new class graph node
                except KeyError:
                    new_label = max(class_graph.nodes()) + 1
                    class_graph.add_node(new_label, nx_graph=new_graph,
                                         edges=new_edges, hash=new_hash)
                    if save_edges:
                        class_graph.add_edge(graph_label, new_label,
                                             equivs=[equiv_nodes],
                                             ops=[op_label])
                    class_graph.member_hash_table.update({new_hash: new_label})
                    queue.append(new_label)
    return class_graph


def int_relabel_graph(graph):
    """
    Relabels graphs with tuple node names to int node names.
    Returns relabelled graph and node mapping applied.
    """
    int_labels = {node: i for i, node in enumerate(graph.nodes())}
    int_graph = nx.relabel_nodes(graph, int_labels)
    return int_graph, int_labels


def explore_lc_orbit(init_graph, save_edges=True, verbose=True):
    """ Explores the LC equivalence class orbit up to isomorphism """
    # Tries to get graph dimensions p^m. If not assigned assumes d = 2
    p = init_graph.__dict__.get('prime', 2)
    m = init_graph.__dict__.get('power', 1)
    if m == 1 and not nx.is_connected(init_graph):
        raise TypeError("Initial graph must be connected.")
    # Creates finite field for arithmetic and maps graph edge int weights
    if m > 1:
        # Creates list of local operations accessible for search
        local_ops = [('CC%d(c,%d)' % (a, t), make_pp_CC_a(a, t))
                     for a in range(1, p) for t in range(m)]
        local_ops += [('EM%d' % (b), make_EM_b(b))
                      for b in range(2, p)]
    elif p > 2 and m == 1:
        local_ops = [('LC' + str(a), make_LC_a(a)) for a in range(1, p)]
        local_ops += [('EM' + str(b), make_EM_b(b)) for b in range(2, p)]
    else:
        local_ops = [('LC', qubit_LC)]
    # Performs orbit search
    class_graph = queued_orbit_search(init_graph, local_ops, save_edges,
                                      verbose)
    # Adds weighted edge data for qudit class graphs
    for node, data in class_graph.nodes.data():
        graph = data['nx_graph']
        edges = data['edges']
        if p ** m > 2:
            data['edges'] = [(u, v, graph[u][v]['weight']) for u, v in edges]
        else:
            data['edges'] = [(u, v, 1) for u, v in edges]
    return class_graph


def get_min_edge_reps(class_graph):
    """ Returns all minimum edge representations for a given LC orbit """
    min_edges = min(len(graph['edges']) for graph in class_graph.node.values())
    min_edge_reps = {key: graph for key, graph in class_graph.node.iteritems()
                     if len(graph['edges']) == min_edges}
    return min_edge_reps


def get_max_edge_reps(class_graph):
    """ Returns all minimum edge representations for a given LC orbit """
    max_edges = max(len(graph['edges']) for graph in class_graph.node.values())
    max_edge_reps = {key: graph for key, graph in class_graph.node.iteritems()
                     if len(graph['edges']) == max_edges}
    return max_edge_reps


def export_class_graph(class_graph, filename, min_edge_reps=False):
    """ Exports class graph to JSON file """
    # Removes all networkx graphs and gets class graph data
    for node, attrs in class_graph.node.iteritems():
        class_graph.node[node].pop('nx_graph')
    class_graph_data = {key: value for key, value
                        in json_graph.node_link_data(class_graph).items()
                        if key in ('nodes', 'links')}
    # Exports class graph to JSON format
    cg_filename = filename + '.json'
    with open(cg_filename, 'w') as fp:
        json.dump(class_graph_data, fp)
    # Finds minimum edge representatives and exports to file
    if min_edge_reps:
        min_edge_reps = get_min_edge_reps(class_graph)
        mer_filename = filename + '_MERs.json'
        with open(mer_filename, 'w') as fp:
            json.dump(min_edge_reps, fp)
    return class_graph_data


def export_class_register(class_graph, filename, min_edge_reps=False):
    """ Exports list of all class members to file """
    # Creates class member register
    register = [[node, attrs['edges'], attrs['hash']]
                for node, attrs in class_graph.node.iteritems()]
    reg_filename = filename + '.csv'
    # Exports register to file
    with open(reg_filename, 'w') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(register)
    # Finds minimum edge representatives and outputs to file
    if min_edge_reps:
        min_edges = len(min(register, key=lambda a: len(a[1]))[1])
        mer_register = [m for m in register if len(m[1]) == min_edges]
        mer_filename = filename + '_MERs.csv'
        with open(mer_filename, 'w') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerows(mer_register)
    return register


if __name__ == '__main__':
#     from psuedo_graphs import gen_psuedo_graph_edge_map, psuedo_to_real, \
#         create_psuedo_graph
#     prime = 2
#     power = 2
#     c_map = gen_psuedo_graph_edge_map(prime, power)
#     pprint(c_map)
#     c_edges = [(0, 1, 1)]
#     psuedo_graph = create_psuedo_graph(c_edges, prime, power, c_map)
#     print psuedo_graph.edges()
#     real_graph = psuedo_to_real(psuedo_graph)
#     print real_graph.edges(data='weight')
#     print real_graph.nodes()
#     class_graph = explore_lc_orbit(real_graph, False)
#     print class_graph.edges()
#     filename = 'class_graphs/tests/ququart_2GHZ_test'
#     register = export_class_register(class_graph, filename, min_edge_reps=True)
#     pprint(register)
#     class_graph_data = export_class_graph(
#         class_graph, filename, min_edge_reps=True)
#     pprint(class_graph_data)

    from graph_builders import create_prime_graph
    prime = 3
    w_edges = [(0, 1, 1), (1, 2, 2)]
    qutrit_g = create_prime_graph(w_edges, prime)
    class_graph = explore_lc_orbit(qutrit_g)
    filename = 'qutrit_ring_class_graph'
    class_graph_data = \
        export_class_graph(class_graph, filename, min_edge_reps=True)
    print class_graph_data
    # pprint(class_graph_data)

    # prime = 3
    # power = 2
    # c_map = gen_psuedo_graph_edge_map(prime, power)
    # # pprint(c_map)
    # c_edges = [(0, 1, 14)]
    # psuedo_graph = create_psuedo_graph(c_edges, prime, power, c_map)
    # real_graph = psuedo_to_real(psuedo_graph)
    # # print real_graph.edges()
    # class_graph = explore_lc_orbit(real_graph)
    # filename = 'class_graphs/tests/qudit_3_2_2GHZ_test'
    # class_graph_data = export_class_graph(
    #     class_graph, filename, min_edge_reps=True)
    # pprint(class_graph_data)

    # prime = 2
    # power = 2
    # c_map = gen_psuedo_graph_edge_map(prime, power)
    # pprint(c_map)
    # c_edges = [(0, 1, 7), (1, 2, 7)]
    # psuedo_graph = create_psuedo_graph(c_edges, prime, power, c_map)
    # real_graph = psuedo_to_real(psuedo_graph)
    # class_graph = explore_lc_orbit(real_graph, save_edges=False)
    # filename = 'class_graphs/tests/ququart_3GHZ_test'
    # class_graph_size = sys.getsizeof(class_graph)
    # print
    # print class_graph_size
    # class_graph_data = export_class_graph(
    #     class_graph, filename, min_edge_reps=True)
    # class_graph_size = sys.getsizeof(class_graph_data)
    # print
    # print class_graph_size
    # pprint(class_graph_data)

    # prime = 2
    # power = 2
    # c_map = gen_psuedo_graph_edge_map(prime, power)
    # # pprint(c_map)
    # c_edges = [(0, 1, 7), (0, 2, 7), (0, 3, 7)]
    # psuedo_graph = create_psuedo_graph(c_edges, prime, power, c_map)
    # real_graph = psuedo_to_real(psuedo_graph)
    # class_graph = explore_lc_orbit(real_graph)
    # filename = 'class_graphs/ququarts/ququart_4GHZ'
    # class_register = export_class_register(class_graph, filename,
    #                                        min_edge_reps=True)
    # class_graph_data = export_class_graph(class_graph, filename,
    #                                       min_edge_reps=True)
    # # pprint(class_graph_data)

    # prime = 3
    # power = 2
    # c_map = gen_psuedo_graph_edge_map(prime, power)
    # # pprint(c_map)
    # c_edges = [(0, 1, 14), (0, 2, 14)]
    # psuedo_graph = create_psuedo_graph(c_edges, prime, power, c_map)
    # real_graph = psuedo_to_real(psuedo_graph)
    # class_graph = explore_lc_orbit(real_graph, save_edges=False)
    # filename = 'class_graphs/qudit_3-2/qudit_3-2_3GHZ'
    # class_register = export_class_register(class_graph, filename,
    #                                        min_edge_reps=True)
    # print len(class_register)
    # class_graph_data = export_class_graph(
    # class_graph, filename, min_edge_reps=True)
    # pprint(class_graph_data)

    # prime = 3
    # power = 2
    # c_map = gen_psuedo_graph_edge_map(prime, power)
    # # pprint(c_map)
    # c_edges = [(0, 1, 14), (0, 2, 14), (0, 3, 14)]
    # psuedo_graph = create_psuedo_graph(c_edges, prime, power, c_map)
    # real_graph = psuedo_to_real(psuedo_graph)
    # class_graph = explore_lc_orbit(real_graph, save_edges=False)
    # filename = 'class_graphs/qudit_3-2/qudit_3-2_4GHZ'
    # class_register = export_class_register(class_graph, filename,
    #                                        min_edge_reps=True)
    # print len(class_register)
    # class_graph_data = export_class_graph(
    #     class_graph, filename, min_edge_reps=True)
    # pprint(class_graph_data)

    # n = 5
    # edges = [(i, (i + 1) % n) for i in range(n)]
    # g = nx.Graph(edges)
    # min_edge_reps = True
    # class_graph = explore_lc_orbit(g)
    # filename = 'class_graphs/tests/5-ring_test'
    # class_graph_data = export_class_graph(
    #     class_graph, filename, min_edge_reps=True)

    # n = 5
    # edges = [(i, (i + 1)) for i in range(n-1)]
    # g = nx.Graph(edges)
    # min_edge_reps = True
    # class_graph = explore_lc_orbit(g)
    # filename = 'class_graphs/tests/5-line_test'
    # class_graph_data = export_class_graph(
    #     class_graph, filename, min_edge_reps=True)

    # n = 6
    # edges = [(i, (i + 1) % n) for i in range(n)]
    # edges += [(0, 2), (5, 3), (1, 4)]
    # print edges
    # g = nx.Graph(edges)
    # min_edge_reps = True
    # class_graph = explore_lc_orbit(g)
    # filename = 'class_graphs/tests/AME6-2_test'
    # class_graph_data = export_class_graph(
    #     class_graph, filename, min_edge_reps=True)

    # from random import randint
    # n = 5
    # prime = 3
    # w_edges = [(i, (i + 1) % n, randint(1, prime-1)) for i in range(n)]
    # g = create_prime_graph(w_edges, prime)
    # min_edge_reps = True
    # class_graph = explore_lc_orbit(g)
    # filename = 'class_graphs/tests/prime_ring%d-%d_test' % (n, prime)
    # class_graph_data = export_class_graph(
    #     class_graph, filename, min_edge_reps=False)

    # from MDS_AME_states import from_MDS_code
    # prime = 2
    # power = 1
    # A = [[1, 1]]
    # graph = from_MDS_code(A, prime, power)
    # class_graph = explore_lc_orbit(graph)
    # filename = 'class_graphs/tests/AME3-2_test'
    # print filename
    # class_graph_data = export_class_graph(
    #     class_graph, filename, min_edge_reps=False)

    # prime = 3
    # power = 1
    # A = [[1, 1], [1, 2]]
    # graph = from_MDS_code(A, prime, power)
    # class_graph = explore_lc_orbit(graph)
    # filename = 'class_graphs/tests/AME4-3_test'
    # print filename
    # class_graph_data = export_class_graph(
    #     class_graph, filename, min_edge_reps=False)

    # prime = 5
    # power = 1
    # A = [[1, 1, 1], [1, 2, 3], [1, 3, 4]]
    # graph = from_MDS_code(A, prime, power)
    # min_edge_reps = True
    # class_graph = explore_lc_orbit(graph)
    # filename = 'class_graphs/tests/AME6-5_test'
    # print filename
    # class_graph_data = export_class_graph(
    #     class_graph, filename, min_edge_reps=False)
