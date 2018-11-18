# Import Python packages
import networkx as nx
# Import local modules
from is_lc_equiv import *

# Create a linear 4 node graph
edges = [(0, 1), (1, 2), (2, 3)]
graph_a = nx.Graph()
graph_a.add_edges_from(edges)

# Create a 4 node ring graph
edges = [(0, 1), (1, 2), (2, 3), (3, 0)]
graph_b = nx.Graph()
graph_b.add_edges_from(edges)

# Create a 4 node ring graph
edges = [(0, 2), (2, 1), (1, 3), (3, 0)]
graph_c = nx.Graph()
graph_c.add_edges_from(edges)

# Checks equivalence between graph A and graph B
is_equiv, local_us = are_lc_equiv(graph_a, graph_b)
print is_equiv, local_us

# Checks equivalence between graph A and graph C
is_equiv, local_us = are_lc_equiv(graph_a, graph_c)
print is_equiv, local_us
