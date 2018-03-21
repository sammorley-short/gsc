# Python packages
import itertools as it
import networkx as nx
import matplotlib.pyplot as plt
from pprint import pprint
# Local modules
from utils import *


def make_crazy(graph, n):
    """ Converts a graph into it's crazy graph encoded version """
    crazy_nodes = {node: [(node, i) for i in range(n)]
                   for node in graph.nodes()}
    crazy_edges = flatten([it.product(crazy_nodes[u], crazy_nodes[v])
                           for u, v in graph.edges()])
    gs = nx.Graph()
    gs.add_nodes_from(flatten(crazy_nodes.values()))
    gs.add_edges_from(crazy_edges)
    if hasattr(g, 'coords'):
        gs.coords = {}
        for node in gs.nodes():
            x, y = coords[node[0]]
            gs.coords.update({node: (x, y - node[1])})
    return gs


if __name__ == '__main__':
    l = 4
    n = 3
    g = nx.path_graph(l)
    coords = {u: (u, 0) for u in g.nodes()}
    g.coords = coords
    pprint(g.edges())
    gs = make_crazy(g, n)
    pprint(gs.edges())
    pprint(gs.coords)
    nx.draw_networkx(gs, pos=gs.coords)
    plt.savefig('graph.png')
