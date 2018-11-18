# Python packages
import itertools as it
import networkx as nx
from pprint import pprint
from random import randint
# Local modules
from utils import *


def random_connected_graph(n):
    """ Generates a Erdos-Renyi G_{n,m} random graph """
    g = nx.Graph([(0, 1), (2, 3)])
    while not nx.is_connected(g):
        edges = randint(n - 1, n * (n - 1) / 2)
        g = nx.dense_gnm_random_graph(n, edges)
    return g


def linear_graph(l):
    """ Creates a linear graph with 2D coordinates """
    g = nx.Graph()
    edges = [((i, 0), (i+1, 0)) for i in range(l-1)]
    g.add_edges_from(edges)
    return g


def make_crazy(graph, n):
    """ Converts a graph into it's crazily encoded version """
    # Defines groupings of crazy nodes and edges between them
    crazy_nodes = {node: [(node, i) for i in range(n)]
                   for node in graph.nodes()}
    crazy_edges = flatten([it.product(crazy_nodes[u], crazy_nodes[v])
                           for u, v in graph.edges()])
    # Creates graph and adds edges and nodes
    crazy_graph = nx.Graph()
    crazy_graph.add_nodes_from(flatten(crazy_nodes.values()))
    crazy_graph.add_edges_from(crazy_edges)
    # Assigns crazy attribute for to_GraphState
    crazy_graph.crazy = True
    return crazy_graph


if __name__ == '__main__':
    l = 5
    g = linear_graph(l)
    cg = make_crazy(g, 6)
    cgs = to_GraphState(cg)
    url = 'https://abv.peteshadbolt.co.uk/SMS_crazy_graphs'
    cgs.url = url
    cgs.push()
