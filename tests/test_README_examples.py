# Import Python packages
import sys
import unittest
import networkx as nx
# Import local modules
from gsc.explore_lc_orbit import explore_lc_orbit, export_class_graph
from gsc.is_lc_equiv import are_lc_equiv
from gsc.graph_builders import create_prime_graph


class TestReadmeExamples(unittest.TestCase):

    def setUp(self):
        pass

    def test_LC_explore_example(self):
        # Create the input graph
        edges = [(0, 1), (1, 2), (2, 3), (3, 4)]
        graph = nx.Graph()
        graph.add_edges_from(edges)
        # Find the class graph
        class_graph = explore_lc_orbit(graph)
        # Export class graph to JSON file
        file_prefix = 'L5_class_graph'
        cg_data = export_class_graph(class_graph, file_prefix)
        # Removes hashes (these may not be the same on a different system)
        for node in cg_data['nodes']:
            node.pop('hash')
        data = {'nodes':
                    [{'edges': [(0, 1, 1), (1, 2, 1), (2, 3, 1), (3, 4, 1)], 'id': 0},
                     {'edges': [(0, 1, 1), (0, 2, 1), (1, 2, 1), (2, 3, 1), (3, 4, 1)], 'id': 1},
                     {'edges': [(0, 1, 1), (1, 2, 1), (1, 3, 1), (2, 3, 1), (3, 4, 1)], 'id': 2},
                     {'edges': [(0, 1, 1), (0, 2, 1), (0, 3, 1), (1, 2, 1), (1, 3, 1), (3, 4, 1)], 'id': 3},
                     {'edges': [(0, 2, 1), (0, 3, 1), (1, 2, 1), (1, 3, 1), (3, 4, 1)], 'id': 4},
                     {'edges': [(0, 2, 1), (0, 3, 1), (0, 4, 1), (1, 2, 1), (1, 3, 1), (1, 4, 1), (3, 4, 1)], 'id': 5},
                     {'edges': [(0, 2, 1), (0, 3, 1), (0, 4, 1), (1, 2, 1), (1, 3, 1), (1, 4, 1), (2, 3, 1), (2, 4, 1)], 'id': 6},
                     {'edges': [(0, 1, 1), (0, 2, 1), (0, 3, 1), (0, 4, 1), (1, 2, 1), (1, 3, 1), (1, 4, 1), (3, 4, 1)], 'id': 7},
                     {'edges': [(0, 1, 1), (0, 2, 1), (0, 3, 1), (0, 4, 1), (2, 3, 1), (2, 4, 1)], 'id': 8},
                     {'edges': [(0, 1, 1), (0, 2, 1), (1, 2, 1), (2, 3, 1), (2, 4, 1), (3, 4, 1)], 'id': 9}],
                'links':
                    [{'source': 0, 'equivs': [[1, 3]], 'target': 1, 'ops': ['LC']},
                     {'source': 0, 'equivs': [[2]], 'target': 2, 'ops': ['LC']},
                     {'source': 1, 'equivs': [[2]], 'target': 8, 'ops': ['LC']},
                     {'source': 1, 'equivs': [[0, 1, 3, 4]], 'target': 9, 'ops': ['LC']},
                     {'source': 2, 'equivs': [[1, 3]], 'target': 3, 'ops': ['LC']},
                     {'source': 3, 'equivs': [[2]], 'target': 4, 'ops': ['LC']},
                     {'source': 3, 'equivs': [[3]], 'target': 5, 'ops': ['LC']},
                     {'source': 4, 'equivs': [[3, 4]], 'target': 8, 'ops': ['LC']},
                     {'source': 4, 'equivs': [[3, 4]], 'target': 7, 'ops': ['LC']},
                     {'source': 5, 'equivs': [[0, 1]], 'target': 6, 'ops': ['LC']},
                     {'source': 5, 'equivs': [[2]], 'target': 7, 'ops': ['LC']},
                     {'source': 6, 'equivs': [[2]], 'target': 9, 'ops': ['LC']},
                     {'source': 7, 'equivs': [[0, 1]], 'target': 8, 'ops': ['LC']}]}
        self.maxDiff = None
        self.assertEqual(data, cg_data)

    def test_LC_equiv_example(self):
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
        self.assertEqual((is_equiv, local_us), (False, None))
        # Checks equivalence between graph A and graph C
        is_equiv, local_us = are_lc_equiv(graph_a, graph_c)
        target_us = [['I', 'H', 'H', 'I'], ['I', 'H', 'SH', 'S'],
                     ['S', 'SH', 'H', 'I'], ['S', 'SH', 'SH', 'S']]
        self.assertEqual((is_equiv, local_us), (True, target_us))

    def test_prime_dimension_explore_example(self):
        # Create the input graph
        prime = 3
        w_edges = [(0, 1, 1), (1, 2, 2)]
        qutrit_g = create_prime_graph(w_edges, prime)
        # Find the class graph
        class_graph = explore_lc_orbit(qutrit_g)
        # Export class graph to JSON file
        filename = 'qutrit_class_graph'
        cg_data = \
            export_class_graph(class_graph, filename)
        # Removes hashes (these may not be the same on a different system)
        for node in cg_data['nodes']:
            node.pop('hash')
        data = {'nodes':
                [{'edges': [(0, 1, 1), (1, 2, 2)], 'id': 0},
                 {'edges': [(0, 1, 2), (1, 2, 2)], 'id': 1},
                 {'edges': [(0, 1, 1), (0, 2, 2), (1, 2, 2)], 'id': 2},
                 {'edges': [(0, 1, 1), (0, 2, 1), (1, 2, 2)], 'id': 3},
                 {'edges': [(0, 1, 1), (1, 2, 1)], 'id': 4},
                 {'edges': [(0, 1, 1), (0, 2, 1), (1, 2, 1)], 'id': 5},
                 {'edges': [(0, 1, 2), (0, 2, 2), (1, 2, 2)], 'id': 6}],
                'links':
                [{'source': 0, 'equivs': [[1]], 'target': 0, 'ops': ['EM2']},
                 {'source': 0, 'equivs': [[0]], 'target': 1, 'ops': ['EM2']},
                 {'source': 0, 'equivs': [[1], [1, 0]], 'target': 2, 'ops': ['LC1', 'LC2']},
                 {'source': 0, 'equivs': [[1], [1, 2]], 'target': 3, 'ops': ['LC2', 'LC1']},
                 {'source': 0, 'equivs': [[2]], 'target': 4, 'ops': ['EM2']},
                 {'source': 1, 'equivs': [[2], [1]], 'target': 2, 'ops': ['LC2', 'LC1']},
                 {'source': 1, 'equivs': [[1]], 'target': 4, 'ops': ['EM2']},
                 {'source': 1, 'equivs': [[1, 0, 2], [1]], 'target': 6, 'ops': ['LC1', 'LC2']},
                 {'source': 2, 'equivs': [[1, 0]], 'target': 2, 'ops': ['EM2']},
                 {'source': 2, 'equivs': [[1, 2], [1, 0]], 'target': 3, 'ops': ['LC2', 'LC1']},
                 {'source': 2, 'equivs': [[1, 0, 2]], 'target': 5, 'ops': ['EM2']},
                 {'source': 2, 'equivs': [[1, 0, 2], [2]], 'target': 6, 'ops': ['LC2', 'LC1']},
                 {'source': 3, 'equivs': [[1, 2]], 'target': 3, 'ops': ['EM2']},
                 {'source': 3, 'equivs': [[1], [0]], 'target': 4, 'ops': ['LC2', 'LC1']},
                 {'source': 3, 'equivs': [[1, 0, 2], [0]], 'target': 5, 'ops': ['LC1', 'LC2']},
                 {'source': 3, 'equivs': [[0]], 'target': 6, 'ops': ['EM2']},
                 {'source': 4, 'equivs': [[1], [1, 0, 2]], 'target': 5, 'ops': ['LC1', 'LC2']}]}


if __name__ == '__main__':
    unittest.main()
