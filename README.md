# gsc (graph-state compass)

Version 2.0

*by Sam Morley-Short*

## Why gsc?

"gsc" or "graph state compass" is a Python package containing a number of tools used for mapping and depicting the local Clifford equivalence classes of quantum graph states.
Specifically, it contains functions that:

* Explores the entire equivalence class, or *orbit* of given graph state (up to isomorphism), of prime and prime-power dimension graph states.
* Detect if two qubit graph states $\left|G\right\rangle$, $\left|G^\prime\right\rangle$ are LC-equivalent and if so returns all local unitaries $U$ such that $\left|G^\prime\right\rangle = U\left|G\right\rangle$.
* Find a graphs' minimum and maximum edge representatives.
* Enumeration of all non-isomorphic equivalence classes for a given $n$-qubit $p^m$-dimensional graph state.
* Visualize graphs and their equivalence classes

Specific instructions on using these tools are given below and a glossary is also provided.

N.B. this README contains latex-formatted equations. If viewing on GitHub these can be easily rendered using browser extensions such as [MathJax for GitHub](https://chrome.google.com/webstore/detail/github-with-mathjax/ioemnmodlmafdkllaclgeombjnmnbima). For mac users, [MacDown](https://macdown.uranusjr.com/) is recommended for offline viewing.

## Installation

To install gsc, clone the repo, and then run `python setup.py install` in the top-level directory.

Note that gsc depends on a forked version of the Pynauty package, which can be found [here](https://github.com/sammorley-short/pynauty).
Pynauty is a Python wrapper for the C program Nauty for computing automorphism groups of graphs.
Please follow the instructions outlined in the linked repo's README for instructions on how to install Pynauty and Nauty.

## Exploring LC-equivalence classes of qubit graph states

### What?

Full exploration of a graph state's local complementation class is not generally efficient.
One reason for this is that within a class there may exist many isomorphs of a single graph, each produced by some unique series of LC's applied to some input graph.
However, in many cases one only cares about the general structure of graphs attainable from the input graph, rather than specific labellings.

By leveraging properties of graphs related to their [automorphisms](https://www.wikiwand.com/en/Graph_automorphism), our algorithm explores the class up to isomorphism faster than previous approaches.
This is achieved using a recursive tree-search approach, whereby each edge that links two graphs within the class (indicating they are separated by at most one LC operation) is only traversed once.
This is acheived by identifying sets of nodes which are equivalent in the graph up to local complementation and only performing LC on one node of each set.

### How?

The search itself is a breadth-first search performed by the function `explore_lc_orbit` found in `explore_class.py`.
This function takes some NetworkX graph as input and outputs a JSON format database of the graph's LC class.

For example, the following snippet of code will find the LC-equivalence class of the 5-qubit linear cluster state:

```python
# Import Python packages
import networkx as nx
# Import local modules
from explore_class import explore_lc_orbit, export_class_graph

# Create the input graph
edges = [(0, 1), (1, 2), (2, 3), (3, 4)]
graph = nx.Graph()
graph.add_edges_from(edges)

# Find the class graph
class_graph = explore_lc_orbit(graph)

# Export class graph to JSON file
filename = 'L5_class_graph.json'
export_class_graph(class_graph, filename)
```

Here the first output that is produced is `class_graph`, a NetworkX graph object whose nodes represent the graphs within the class and edges the local complementations that connect them.
This object is then exported to JSON using a standard NetworkX format that can be reloaded into NetworkX to reproduce the class graph later.
Within the JSON dictionary there are two important data structures representing the class graph's nodes and edges respectively.
Firstly, the list keyed by `"nodes"` stores a list of dictionaries, each representing a node (i.e. graph) in the orbit containing the following information:

* `"edges"`: The edges of the graph the node represents.
* `"hash"`: A hash value of the canonically relabelled version of the graph (such that two isomorphic graphs will have the same hash value).
* `"id"`: A short integer label for the graph used to define the class graph's edges.

Secondly, the list keyed by `"links"` gives a list of dictionaries representing the edges of the class graph.
Each dictionary representing the edge $(i, j)$ contains a `"source": i` and `"target": j` as well as operations were applied to which nodes to map the source to target graph state, keyed by `"equivs"` and `"ops"` respectively.
For qubit graphs `"equivs"` and `"ops"` are single element lists because only a single operation exists---namely, local complementation---however this is not the case for higher-dimensional graph states.

**Notes:**

* All node labels on the input graph must be integers. If this is not the case, `int_relabel_graph` can be used to return a relabelled graph and the labelling applied.
* Local complementation operations on degree one qubits are trivial operations and so are ignored.
* `explore_lc_class` has the following optional keyword arguments:
	* `save_edges=True`: if set to `False` then only the class members themselves are stored in the class graph, discarding the operations that connect them.
		For example, this is used for member enumeration, where the class' specific structure is not sought.
	* `verbose=True`: By default the ratio of explored graphs to known members is displayed during the search.
		To turn this off set `verbose=False`.
* `export_class_graph` has the optional keyword argument `min_edge_reps=False`. Setting `min_edge_reps=True` when also export a list of the classes Minimum Edge Representatives (MERs).
* `export_class_register` can be used to export a register of class members to a CSV table containing each graph's ID, edge list and hash.
	This should be used when the equivalence class' internal structure is not needed, such as in enumeration.

## Testing for LC-equivalence

### What?

For the reduced task of testing whether two known graphs $\left|G\right\rangle$ and $\left|G^\prime\right\rangle$ are LC-equivalent, we include an implementation of an algorithm [originally presented](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.70.034302) by Van Den Nest, et. al.
In the case that two graphs are equivalent, the algorithm also returns all local unitaries $U$ such that $\left|G^\prime\right\rangle = U\left|G\right\rangle$.

### How?

LC equivalence is checked by function `are_lc_equiv` found in `is_lc_equiv.py`. The function takes as input two NetworkX graphs and outputs the tuple `(are_lc_equiv, local_ops)` which are a boolean and a list of valid unitaries.

For example, consider checking the local equivalence between some set of 4-qubit graphs:

```python
# Import Python packages
import networkx as nx
# Import local modules
from is_lc_equiv import are_lc_equiv

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

## Higher-dimension graph states

`gsc` can also explore the equivalence classes of higher-dimensional graph states.

### Prime dimension

Quantum stabilizer states of prime local dimension can also be described within the graph-state formalism, as described in the following [paper](https://arxiv.org/abs/quant-ph/0610267).
As such, by including the higher-dimensional local complementation and edge multiplication operations, `gsc` can also explore the equivalence classes of prime dimension graph states.

For example, the following snippet can be used to export the equivalence class of a three-qutrit state:

```python
# Import the prime graph state builder and class explorer and exporter
from graph_builders import create_prime_graph
from explore_lc_orbit import explore_lc_orbit, export_class_graph
# Create the input graph
prime = 3
w_edges = [(0, 1, 1), (1, 2, 2)]
qutrit_g = create_prime_graph(w_edges, prime)
# Find the class graph
class_graph = explore_lc_orbit(qutrit_g)
# Export class graph to JSON file
filename = 'qutrit_class_graph'
class_graph_data = export_class_graph(class_graph, filename)
```
where `class_graph_data` contains the exported JSON-formatted dictionary.

### Prime-power dimension

Quantum stabilizer states of prime-power dimension can also be mapped onto prime-dimensional graph states and therefore also described within the graph-state formalism.
Under the isomorphism described in Chapter 5 of my PhD thesis, each qudit of local dimension $p^m$ (for $p$ prime and some integer $m>2$) is mapped to a party of $m$ $p$-dimensional qudits between which entangling operations are now considered local.
In this case, the local complementation operation is replaced with the so-called controlled complementation operation, which with edge multiplcation completes the set of local graph-state operations.

```python
# Import the pseudo graph state builder and class explorer and exporter
from psuedo_graphs import gen_psuedo_graph_edge_map, create_psuedo_graph, psuedo_to_real
from explore_lc_orbit import explore_lc_orbit, export_class_graph
# Create the initial prime-power pseudo graph state
prime, power = 2, 2
c_map = gen_psuedo_graph_edge_map(prime, power)
c_edges = [(0, 1, 5), (1, 2, 2)]
psuedo_graph = create_psuedo_graph(c_edges, prime, power, c_map)
# Convert prime-power pseudo graph state to "real" prime graph state
real_graph = psuedo_to_real(psuedo_graph)
# Find the class graph
class_graph = explore_lc_orbit(real_graph, save_edges=False)
# Export class graph to JSON file
filename = 'p2_m2_qudit_class_graph'
class_graph_data = export_class_graph(class_graph, filename)
# class_graph_data also contains the exported JSON-formatted dictionary
```
where `class_graph_data` contains the exported JSON-formatted dictionary.

Since prime-power graphs often contain many, many edges, for convenience, we have initialised the prime-power graph state using its *pseudo graph state* representation.
A pseudo graph state is a representation of the prime-power graph state using a directed graph where each edge label or *colour* represents some configuration of edges between two qubit families on the prime-dimensional graph as defined by some coloured edge mapping.
In the above example, this mapping is defined by:

```python
c_map = {0: [],
         1: [(0, 0, 1)],
         2: [(0, 1, 1)],
         3: [(1, 0, 1)],
         4: [(1, 1, 1)],
         5: [(0, 0, 1), (0, 1, 1)],
         6: [(0, 0, 1), (1, 0, 1)],
         7: [(0, 0, 1), (1, 1, 1)],
         8: [(0, 1, 1), (1, 0, 1)],
         9: [(0, 1, 1), (1, 1, 1)],
         10: [(1, 0, 1), (1, 1, 1)],
         11: [(0, 0, 1), (0, 1, 1), (1, 0, 1)],
         12: [(0, 0, 1), (0, 1, 1), (1, 1, 1)],
         13: [(0, 0, 1), (1, 0, 1), (1, 1, 1)],
         14: [(0, 1, 1), (1, 0, 1), (1, 1, 1)],
         15: [(0, 0, 1), (0, 1, 1), (1, 0, 1), (1, 1, 1)]}
```
where each edge $(i, j, w)$ defines the edge between family members $i$ and $j$ of the source and sink nodes of weight $w$.
For example, in the above example the pseudo edge list `[(0, 1, 5), (1, 2, 2)]` represents the prime-dimension graph with edges `[((0, 0), (1, 0)), ((0, 0), (1, 1)), ((1, 0), (2, 1))]` (and nodes `[(0, 0), (0, 1), (1, 0), (1, 1), (2, 0), (2, 1)]`), where unit edge weights have been omitted for brevity.
It is easy to see that while such states are convenient to reduce excessive listing of edges, they are clearly not unique.
For example, the above state could have been equally represented by the pseudo edge list `[(1, 0, 6), (2, 1, 3)]` (recall that the graph is directed, and so the edge $(i, j, c)$ represents the edge $i \rightarrow j$ of colour $c$).

While the pseudo graph state representation is purely used for convenience, it becomes useful for class enumeration, where listing all the edges of large states becomes unweildy.

## Class enumeration

`gsc` also performs equivalence class enumeration, and is performed as follows.
Firstly, initialise a register of all possible $p^m$-dimension $n$-qudit graph states (in `remaining_graphs.csv` with parameters $p, m, n$ stored in `state_params.csv`), an empty hash lookup table of found graphs (in `graph_hashes.csv`) and an empty class database directiory (in `/classes`).
In the register, a graph is stored as a list or "configuration" of pseudo-edges indexed by list of edges stored in `edge_index.csv` and are associated to their associated prime-dimensional graph-state via the pseudo-edge map defined in `psuedo_edge_map.csv`.
Next, for each state:

1. **If disconnected:** Remove graph and all isomorphs from graph register.
2. **If graph hash already in hash table:** Remove graph and all isomorphs from graph register.
3. **Else:** **i)** Explore local equivalence class to find all members, **ii)** Export class to database, **iii)** Find all isomorphs of all members and their hashes, **iv)** Remove all isomorphs from register and add all hashes to hash table.

This is repeated until the graph register is empty.

**Notes:**

* Each output equivalence class is named by the edge configuration of their first member and contains a register of each graph it contains.
	Each graph in a class register is stored as a row of three columns, denoting the graph's ID, prime-dimensional edge list, and it's hash value.
	For $p=2$ graph states, edge weights are omitted, and for $m=1$ graph states, so are family member labels.
* While it may appear that finding isomorphs may be redundant in step 2), given that they are also found in step 3).
	This is because graphs states are stored in their pseudo-edge representation (which saving on reading and writing entire edge lists associated with their prime-dimensional representation).
	Hence, some graphs with differing pseudo-edge representation actually represent the same state in their prime-dimension.
	For example, the 3-ququart states with edge lists:

	```
	[((0, 0), (1, 0)), ((1, 0), (2, 0)), ((2, 0), (0, 0))] and
	[((0, 1), (1, 1)), ((1, 1), (2, 1)), ((2, 1), (0, 1))]
	```

	where unit edge weights have been omitted for brevity).
	Clearly both edge lists represent equivalent states, yet they have weighted pseudo edge lists of

	```
	[(0, 1, 1), (1, 2, 1), (2, 0, 1)]
	[(0, 1, 4), (1, 2, 4), (2, 0, 4)]
	```
	respectively.
	Because isomorphic configurations are found while in the pseudo-edge representation, the following two states will not be distinguished until their prime-dimensional graph states are hashed.

## Dependancies

This module relies on the following packages:

* Nauty (found [here](http://pallini.di.uniroma1.it/))
* Pynauty-hack (found [here](https://github.com/sammorley-short/pynauty-hack))
* Various python modules: NetworkX, Numpy, tqdm, abp, matplotlib (which should all installed via pip at installation)

## Glossary

* **Graph state**: a quantum state $\left|G\right\rangle$ represented by some graph $G$.
* **Local complementation (LC)**: the graph operation that represents taking graph states to other graph states using only local Clifford unitaries.
* **Class graph**: a simple, connected and undirected graph representing the structure of an LC-equivalence class. Each node represents some graph state within the class with each edge between two nodes the application of an LC operation that takes one graph to another.
* **Isomorphic**: Two graphs are isomorphic if they are the same up to relabelling. This is equivalent to saying they have the same [*canonical labelling*](https://www.wikiwand.com/en/Graph_canonization).
* **Automorphism**: An automorphism of a graph is any relabelling that produces the same graph. E.g. for a ring graph with $V = \\{0, \ldots, n\\}$ nodes, the node relabelling $i \mapsto (i+1) \bmod{|V|}$ is a valid automorphism as it produces the same graph.
* **NetworkX**: A useful python package for creating and manipulating graphs. On Unix systems it can be most easily installed via pip through the command `$ pip install networkx`. Documentation can be found [here](https://networkx.github.io/documentation/stable/index.html).
* **ABP**: A python package for efficiently simulating and visualising quantum graph states based on Anders and Briegel's original algorithm. It can also be installed via pip, with installation instructions and documentation can be found [here](https://github.com/peteshadbolt/abp).
* **Nauty**: A popular C library for finding graph auto- and isomorphisms, found [here](http://pallini.di.uniroma1.it/) which must be installed. As well as this, our code relies on a hacked version of the python wrapper `pynauty` found [here](https://github.com/sammorley-short/pynauty-hack).