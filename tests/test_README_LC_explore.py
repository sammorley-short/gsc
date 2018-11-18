# Import Python packages
import networkx as nx
# Import local modules
from explore_class import *

# Create the input graph
edges = [(0, 1), (1, 2), (2, 3), (3, 4)]
graph = nx.Graph()
graph.add_edges_from(edges)

# Find the class graph
class_graph = explore_LC_isomorphic_orbit(graph)

# Export class graph to JSON file
filename = 'L5_class_graph.json'
export_class_graph(class_graph, filename)
