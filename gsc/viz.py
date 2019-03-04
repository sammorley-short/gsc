# Local modules
from gsc.utils import to_GraphState


def push_graph_to_abv(graph, url=None):
    """
        Converts a graph to a ABP GraphState and pushes to abv server.
        Returns GraphState for debugging
    """
    gs = to_GraphState(graph)
    if url is not None:
        gs.url = url
    gs.push()
    return gs
