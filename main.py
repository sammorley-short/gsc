# Local modules
from viz import *
from explore_class import *
from graph_builders import *


if __name__ == '__main__':
    url = "https://abv.peteshadbolt.co.uk/SMS_crazy_graph-main"
    n = 4
    l = 5
    lin_cg = make_crazy(linear_graph(l), n)
    lin_cgs = push_graph_to_abv(lin_cg, url)
    edge = (((2, 0), 0), ((1, 0), 0))
    lin_cg = edge_lc(lin_cg, edge)
    lin_cgs = push_graph_to_abv(lin_cg, url)

    # filename = 'class_graphs/crazy_graphs/lin_L%d_N%d.json' % (l, n)
    # relab_lin_cg, int_relab = int_relabel_graph(lin_cg)
    # class_graph = explore_LC_isomorphic_orbit(relab_lin_cg)
    # class_graph_data = export_class_graph(class_graph, filename)
    # pprint(class_graph_data)
