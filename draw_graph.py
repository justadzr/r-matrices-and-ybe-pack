import sympy as sp, triple, ybe, multiprocessing as mp
import networkx as nx, matplotlib.pyplot as plt

def red(a, b):
    return (a - 1) % b + 1

n=12
with open(f"nonassociative-affine-triples-{n}.txt", "r") as f:
    unclean_triples = f.read()[2:-2].split('], [')
triples = []
for unclean in unclean_triples:
    t = str.maketrans({',': ' ', '[': ' ', ']': ' '})
    tuple = list(map(int, unclean.translate(t).split()))
    if tuple:
        triples += [triple.BDTriple(tuple)]

trip = triple.BDTriple(tuple)

counter = 0
for trip in triples:
    counter += 1
    tup = trip.tuple
    connected_components = trip.connected_components()
    connected_components_img = trip.connected_components_img()
    edges = []
    G = nx.DiGraph()
    for i in range(n):
        G.add_edge(f"{i+1}", f"{red(i+2, n)}", color='w', weight=0)
    node_color = [-1] * n
    num_comp = len(connected_components)

    for i in range(num_comp):
        component = connected_components[i]
        component_img = connected_components_img[i]
        for j in range(len(component)):
            alpha = component[j]
            beta = trip.T(alpha)
            G.add_edge(f"{alpha}", f"{beta}", color='b', weight=2, node_color=i+1)
            node_color[alpha-1] = i

    pos = nx.circular_layout(G)
    colors = nx.get_edge_attributes(G,'color').values()
    weights = nx.get_edge_attributes(G,'weight').values()
    nx.draw(G, pos=pos, with_labels=True, edge_color=colors, width=list(weights), node_color=node_color, cmap=plt.cm.tab10)
    plt.savefig(f'graph{n}/{counter}-graph.png')
    plt.close()