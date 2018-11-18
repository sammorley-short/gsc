# Tools for exploring graph state equivalence classes

*By Sam Morley-Short*

**Version 1.0** 

---

This repository contains a number of tools used for exploring the local complementation equivalence classes of quantum graph states.
Specifically, it contains functions that:

* Explores the entire LC-equivalence class, or *orbit* of a given graph state (up to isomorphism).
* Detect if two graph states $\left|G\right\rangle$, $\left|G^\prime\right\rangle$ are LC-equivalent and if so returns all local unitaries $U$ such that $\left|G^\prime\right\rangle = U\left|G\right\rangle$.
* Find a graphs' minimum and maximum edge representatives.
* Visualize graphs and their equivalence classes

Specific instructions on using these tools are given below and a glossary is also provided.

N.B. this README contains latex-formatted equations. If viewing on GitHub these can be easily rendered using browser extensions such as [MathJax for GitHub](https://chrome.google.com/webstore/detail/github-with-mathjax/ioemnmodlmafdkllaclgeombjnmnbima). For mac users, [MacDown](https://macdown.uranusjr.com/) is recommended for offline viewing. 

## Exploring LC-equivalence classes

### What is?

Full exploration of a graph's local complementation class is not generally efficient.
One reason for this is that within a class there may exist many isomorphs of a single graph, each produced by some unique series of LC's applied to some input graph.
However, in many cases one only cares about the general structure of graphs attainable from the input graph, rather than specific labellings.

By leveraging properties of graphs related to their [automorphisms](https://www.wikiwand.com/en/Graph_automorphism), our algorithm efficiently explores the class up to isomorphism.
This is achieved using a recursive tree-search approach, whereby each edge that links two graphs within the class (indicating they are separated by at most one LC operation) is only traversed once, and hence the class graph 

### How do?

The search itself is performed by the function `explore_LC_isomorphic_orbit` found in `explore_class.py`.
This function takes some NetworkX graph as input and outputs a JSON format database of the graph's LC class.

In practise, this is achieved by the recipe depicted in the following python script (found in `tests/test_README_LC_explore`:

```python
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
```

For example, the above code finds the orbit for the 5-qubit linear graph state.
The output JSON file uses a standard NetworkX format so that in principle can be reloaded into NetworkX to reproduce the class graph in some other code.
However, from our perspective it contains two important data structures representing the class graph's nodes and edges respectively.
Firstly, the list keyed by `"nodes"` stores a list of dictionaries, each representing a node (i.e. graph) in the orbit containing the following information:

* `"edges"`: The edges of the graph the node represents.
* `"hash"`: A hash value of the canonically relabelled version of the graph (such that two isomorphic graphs will have the same hash value).
* `"id"`: A short integer label for the graph used to define the edges.

Secondly, the list keyed by `"adjacency"` gives an list in which the $i$th element is a list of dictionaries that represent each edge $(i, j)$ connected to node $i$.
Each dictionary representing $(i, j)$ contains the `"id"` of node $j$ and `"equiv"`, a list of nodes in the graph $i$ that produce the same graph $j$ after LC (up to isomorphism).

**Notes:**

* All node labels on the input graph must be integers. If this is not the case, `int_relabel_graph` can be used to return a relabelled graph and the labelling applied. 

## Testing for LC-equivalence

### What is?

For the reduced task of testing whether two known graphs $\left|G\right\rangle$ and $\left|G^\prime\right\rangle$ are LC-equivalent, we include an implementation of an algorithm [originally presented](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.70.034302) by Van Den Nest, et. al.
In the case that two graphs are equivalent, the algorithm also returns all local unitaries $U$ such that $\left|G^\prime\right\rangle = U\left|G\right\rangle$.

### How do?

LC equivalence is checked by function `are_lc_equiv` found in `is_lc_equiv.py`. The function takes as input two NetworkX graphs and outputs the tuple `(are_lc_equiv, local_ops)` which are a boolean and list of possible $U$.

For example, consider the following python script (found in `tests/test_README_LC_equiv.py`):

```python
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
```

which outputs the following:

```
False None
True [['I', 'H', 'H', 'I'], ['I', 'H', 'SH', 'S'], ['S', 'SH', 'H', 'I'], ['S', 'SH', 'SH', 'S']]
```

The validity of the output can be seen by referring to Fig. 7 of the following [paper](https://arxiv.org/abs/quant-ph/0307130v7).

## Dependancies

This module relies on the following packages:

* Nauty (found [here](http://pallini.di.uniroma1.it/))
* Pynauty-hack (found [here](https://github.com/sammorley-short/pynauty-hack))
* Various python modules: NetworkX, Numpy, pprint, abp (all installed via pip)

## Glossary

* **Graph state**: a quantum state $\left|G\right\rangle$ represented by some graph $G$.
* **Local complementation (LC)**: the graph operation that represents taking graph states to other graph states using only local Clifford unitaries.
* **Class graph**: a simple, connected and undirected graph representing the structure of an LC-equivalence class. Each node represents some graph state within the class with each edge between two nodes the application of an LC operation that takes one graph to another.
* **Isomorphic**: Two graphs are isomorphic if they are the same up to relabelling. This is equivalent to saying they have the same [*canonical labelling*](https://www.wikiwand.com/en/Graph_canonization).
* **Automorphism**: An automorphism of a graph is any relabelling that produces the same graph. E.g. for a ring graph with $V = \{0, \ldots, n\}$ nodes, the node relabelling $i \mapsto (i+1) \;\textrm{mod} \; |V|$ is a valid automorphism as it produces the same graph.
* **NetworkX**: A useful python package for creating and manipulating graphs. On Unix systems it can be most easily installed via pip through the command `$ pip install networkx`. Documentation can be found [here](https://networkx.github.io/documentation/stable/index.html).
* **ABP**: A python package for efficiently simulating and visualising quantum graph states based on Anders and Briegel's original algorithm. It can also be installed via pip, with installation instructions and documentation can be found [here](https://github.com/peteshadbolt/abp).
* **Nauty**: A popular C library for finding graph auto- and isomorphisms, found [here](http://pallini.di.uniroma1.it/) which must be installed. As well as this, our code relies on a hacked version of the python wrapper `pynauty` found [here](https://github.com/sammorley-short/pynauty-hack).