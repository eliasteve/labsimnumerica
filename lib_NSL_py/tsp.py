import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

def plot_network(G, title):
    """
    Plots the network.
    
    Args:
        G: Graph to plot
        title: title for the plot
    """
    pos = nx.get_node_attributes(G, 'pos')
    nx.draw(G, pos, with_labels=False, node_size=10, node_color='skyblue', font_size=10, font_weight='bold', arrows=True)
    plt.title(title)
    plt.show()



def create_graph(cities, path):
    """
    Given a list of cities and a path through them, creates a graph containing the path.

    Args:
        cities: list of cities
        path: path through the cities as a list of city labels (here numbers)
    """
    G = nx.DiGraph()
    for n, pos in enumerate(cities):
        G.add_node(n, pos=pos)
    
    for n1, n2 in zip(path, path[1:]):
        G.add_edge(n1, n2)

    return G



def compute_cycle_length(G, cycle):
    """
    Computes the length of the cycle in the graph.
    
    Args:
        G: graph
        cycle: cycle whose length is to be computed
    """
    length = 0
    for (id1, id2) in zip(cycle, cycle[1:]):
        pos1 = G.nodes[id1]['pos']
        pos2 = G.nodes[id2]['pos']
        distance = np.sqrt((pos2[0] - pos1[0])**2 + (pos2[1] - pos1[1])**2)
        length += distance
    return length
