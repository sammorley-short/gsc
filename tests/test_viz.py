# Python packages
import sys
import ast
import unittest
import networkx as nx
from abp import GraphState
# Local modules
sys.path.append('..')
import viz
from graph_builders import linear_graph, make_crazy
from utils import to_GraphState

url = "https://abv.peteshadbolt.co.uk/beam-robot-reckless-facade"


def encode_dict(d, codec='utf8'):
    """ Processes a dictionary loaded from JSON """
    ks = d.keys()
    for k in ks:
        val = d.pop(k)
        if isinstance(val, unicode):
            val = val.encode(codec)
        elif isinstance(val, dict):
            val = encode_dict(val, codec)
        if isinstance(k, unicode):
            k = k.encode(codec)
        d[k] = val
    return d


def pull_GraphState_nodes_edges(url):
    """
    Pulls GraphState from url and processes nodes & edges to standard form
    """
    gs2 = GraphState()
    gs2.pull(url)
    gs2_nodes = encode_dict(gs2.node)
    gs2_edges = sorted([(ast.literal_eval(u), ast.literal_eval(v))
                        for u, v in gs2.edgelist()])
    return gs2_nodes, gs2_edges


class TestGraphViz(unittest.TestCase):
    """ Tests pushing to ABV for NetworkX graph vizualisation """

    def setUp(self):
        pass

    def test_linear_graph_push_to_abv(self):
        """ Tests pushing a linear graph to abv and then pulling it back """
        # Creates linear graph
        g = linear_graph(10)
        # Pushes to abv
        viz.push_graph_to_abv(g, url=url)
        # Checks pulled graph same as pushed
        gs1 = to_GraphState(g)
        gs1_nodes = {str(key): value for key, value in gs1.node.items()}
        gs1_edges = sorted(gs1.edgelist())
        # Gets pushed graph state nodes & edges
        gs2_nodes, gs2_edges = pull_GraphState_nodes_edges(url)
        self.assertEqual(gs1_nodes, gs2_nodes)
        self.assertEqual(gs1_edges, gs2_edges)

    def test_crazy_linear_graph_push_to_abv(self):
        """
        Tests pushing a crazy linear graph to abv and then pulling it back
        """
        g = linear_graph(10)
        cg = make_crazy(g, 10)
        # Pushes to abv
        viz.push_graph_to_abv(cg, url=url)
        # Checks pulled graph same as pushed
        gs1 = to_GraphState(cg)
        gs1_nodes = {str(key): value for key, value in gs1.node.items()}
        gs1_edges = sorted(gs1.edgelist())
        # Gets pushed graph state nodes & edges
        gs2_nodes, gs2_edges = pull_GraphState_nodes_edges(url)
        self.assertEqual(gs1_nodes, gs2_nodes)
        self.assertEqual(gs1_edges, gs2_edges)


if __name__ == '__main__':
    unittest.main()
