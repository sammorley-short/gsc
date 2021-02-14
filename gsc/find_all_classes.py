# Python modules
import os
import csv
import sys
import numpy as np
import networkx as nx
import itertools as it
from tqdm import tqdm
from pprint import pprint
# Local modules
from gsc.get_nauty import hash_graph
from gsc.psuedo_graphs import (
    real_graph_to_psu_edges,
    gen_psuedo_graph_edge_map,
    create_psuedo_graph,
    psuedo_to_real,
)
from gsc.explore_lc_orbit import explore_lc_orbit


def init_search_database(prime, power, nodes):
    """ Initialises database for class search """
    # Initialises database folders
    directory = 'class_databases/' + \
        'prime_power_p%d_m%d_n%d' % (prime, power, nodes)
    if not os.path.exists(directory):
        os.makedirs(directory)
    if not os.path.exists(directory + '/classes'):
        os.makedirs(directory + '/classes')
    # Writes edge index (used for edge configurations) to file
    all_edges = list(it.combinations(range(nodes), 2))
    filename = directory + '/edge_index.csv'
    with open(filename, 'w') as file:
        writer = csv.writer(file)
        writer.writerows(all_edges)
    # Writes the state's parameters to file
    filename = directory + '/state_params.csv'
    with open(filename, 'w') as file:
        writer = csv.writer(file)
        writer.writerow([prime, power, nodes])
    # Writes complete list of edge configs to file
    max_edges = nodes * (nodes - 1) // 2
    all_edge_configs = it.product(range(prime ** (power ** 2)),
                                  repeat=max_edges)
    edge_configs = \
        it.ifilter(lambda a: (max_edges - a.count(0)) >= (nodes - 1),
                   all_edge_configs)
    edge_configs = np.fromiter(it.chain.from_iterable(edge_configs), 'int8')
    edge_configs = edge_configs.reshape(-1, max_edges)
    filename = directory + '/remaining_graphs.csv'
    np.savetxt(filename, edge_configs, fmt='%d', delimiter=',')
    # If prime-power, writes psuedo-edge map to file
    if power > 1:
        c_map = gen_psuedo_graph_edge_map(prime, power).items()
        filename = directory + '/psuedo_edge_map.csv'
        with open(filename, 'w') as file:
            writer = csv.writer(file)
            writer.writerows(c_map)
    # Creates graph hash directory
    open(directory + '/graph_hashes.csv', 'w').close()
    return directory


def get_next_graph(directory):
    """ Pops the next graph to analyse from top of remaining_graphs.txt """
    filename = directory + '/remaining_graphs.csv'
    # Reads first line and writes rest to temp. file
    with open(filename) as file:
        edge_config = file.readline()
    # Formats edge config to tuple
    if edge_config:
        edge_config = [int(i) for i in edge_config[:-1].split(',')]
    return edge_config


def remove_found_graphs(directory, edge_configs):
    """ Removes any found graphs from remaining_graphs.txt """
    edge_configs = [map(str, config) for config in edge_configs]
    filename = directory + '/remaining_graphs.csv'
    tmp_filename = directory + '/remaining_graphs_tmp.csv'
    # Reads configs from file and writes any not in found configs to temp.
    with open(filename) as file, open(tmp_filename, "w") as tmp_file:
        reader = csv.reader(file)
        writer = csv.writer(tmp_file)
        for config in reader:
            if config in edge_configs:
                continue
            writer.writerow(config)
    # Removes original file and renames temp. to original
    os.remove(filename)
    os.rename(tmp_filename, filename)


def write_hashes(directory, hashes):
    """
    Writes a sequence of graph hashes to graph_hashes.csv if they haven't
    already been found.
    """
    filename = directory + '/graph_hashes.csv'
    # Reads in found hashes
    with open(filename, 'r') as file:
        reader = csv.reader(file)
        found_hashes = set([int(h[0]) for h in reader])
    # If none have been found, write hashes to file
    hashes = set(hashes)
    if hashes - found_hashes == hashes:
        hashes = map(lambda h: str(h) + '\n', hashes)
        with open(filename, 'a+') as file:
            file.writelines(hashes)
    # If only a subset of the hashes have been found raise an error
    elif hashes - found_hashes != set([]):
        # pprint(hashes)
        pprint(hashes - found_hashes)
        pprint(hashes - (hashes - found_hashes))
        raise Exception("Error: Only some of hashes already known")
        sys.exit()


def found_hash(directory, target_hash):
    """ Returns whether hash has already been found """
    filename = directory + '/graph_hashes.csv'
    with open(filename, 'r') as file:
        reader = csv.reader(file)
        hashes = [int(h[0]) for h in reader]
    return target_hash in hashes


def make_isomorph_func(edge_index, n):
    """ Returns function that makes isomorphic graph edge configurations """
    node_maps = [{i: j for i, j in enumerate(perm)}
                 for perm in it.permutations(range(n), n)]
    iso_indices = [[tuple(sorted([n_map[i], n_map[j]])) for i, j in edge_index]
                   for n_map in node_maps]
    config_perms = [[edge_index.index(edge) for edge in index]
                    for index in iso_indices]

    def isomorph_configs(edge_config):
        iso_configs = [tuple(edge_config[i] for i in perm)
                       for perm in config_perms]
        iso_configs = map(list, set(iso_configs))
        return iso_configs

    return isomorph_configs


def remove_disconnected_configs(directory, edge_config, isomorph_configs):
    """ Removes all graphs which are siliarly disconnected incl. isomorphs """
    # Finds the edge occupancies of all isomorphic configurations
    iso_occs = set(tuple(map(bool, config))
                   for config in isomorph_configs(edge_config))
    iso_occs = map(list, list(iso_occs))
    filename = directory + '/remaining_graphs.csv'
    tmp_filename = directory + '/remaining_graphs_tmp.csv'
    # Reads configs from file and writes any not in found configs to temp.
    with open(filename) as file, open(tmp_filename, "w") as tmp_file:
        reader = csv.reader(file)
        writer = csv.writer(tmp_file)
        for config in reader:
            config_occ = [bool(int(i)) for i in config]
            if config_occ in iso_occs:
                continue
            writer.writerow(config)
    # Removes original file and renames temp. to original
    os.remove(filename)
    os.rename(tmp_filename, filename)


def find_all_classes(directory, power, prime):
    """ Finds all members of all classes """
    # Gets edge indices and state params and generates edge map
    filename = directory + '/edge_index.csv'
    with open(filename, 'r') as file:
        reader = csv.reader(file)
        edge_index = [tuple(map(int, edge)) for edge in reader]
    filename = directory + '/state_params.csv'
    with open(filename, 'r') as file:
        reader = csv.reader(file)
        p, m, n = map(int, next(reader))
    c_map = gen_psuedo_graph_edge_map(p, m)
    # Creates function to produce config isomorphs
    isomorph_configs = make_isomorph_func(edge_index, n)
    pprint(c_map)
    pprint(edge_index)
    # Initialises progress bar
    rg_file = directory + '/remaining_graphs.csv'
    rem_graphs_size = os.path.getsize(rg_file)
    pbar = tqdm(total=rem_graphs_size)
    while True:
        rem_update = rem_graphs_size - os.path.getsize(rg_file)
        if rem_update >= 0:
            pbar.update(rem_update)
        rem_graphs_size = os.path.getsize(rg_file)
        # Waits until a graph is available to process
        edge_config = get_next_graph(directory)
        if not edge_config:
            break
        tqdm.write("Psuedo edge config: %s" % (edge_config,))
        # Create initial graph
        c_edges = [(u, v, w) for (u, v), w in zip(edge_index, edge_config)]
        init_graph = create_psuedo_graph(c_edges, p, m, c_map)
        # Checks if graph is connected
        if not nx.is_connected(init_graph.to_undirected()):
            tqdm.write("Disconnected. Removing isomorphs... ")
            # Removes any isomorphic graphs from remaining
            remove_disconnected_configs(directory, edge_config, isomorph_configs)
            tqdm.write("Done")
            continue
        init_graph = psuedo_to_real(init_graph)
        # Check if graph has already been found for hash test
        graph_hash = hash_graph(init_graph)
        if found_hash(directory, graph_hash):
            tqdm.write("Already found %d" % graph_hash)
            iso_configs = isomorph_configs(edge_config)
            remove_found_graphs(directory, iso_configs)
            continue
        # Explore class graph
        tqdm.write("Exploring class...")
        class_graph = explore_lc_orbit(init_graph, False, False)
        class_register = [[node, attrs['edges'],
                           attrs['hash'], attrs['nx_graph']]
                          for node, attrs in class_graph.node.items()]
        nodes, edges, hashes, graphs = zip(*class_register)
        # Formats edge list based on state parameters
        if power == 1:
            edges = [[(u, v, c) for (u, i), (v, j), c in edge_set]
                     for edge_set in edges]
        if prime == 2:
            edges = [[(u, v) for u, v, c in edge_set]
                     for edge_set in edges]
        # Exports class_register to file
        tqdm.write("Writing register to file...")
        class_register = zip(nodes, edges, hashes)
        config_label = '_'.join(map(str, edge_config))
        filename = directory + '/classes/' + config_label + '.csv'
        with open(filename, 'w') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerows(class_register)
        # Finds and removes any isomorphs
        tqdm.write("Removing isomorphs...")
        psu_edges = [real_graph_to_psu_edges(graph, c_map, edge_index)
                     for graph in graphs]
        edge_configs = [[c for u, v, c in w_edges] for w_edges in psu_edges]
        iso_configs = [iso_config for config in edge_configs
                       for iso_config in isomorph_configs(config)]
        remove_found_graphs(directory, iso_configs)
        # Adds hashes to found hash directory
        write_hashes(directory, hashes)
        tqdm.write("Done")
    pbar.close()
