import sympy as sp, triple, ybe, multiprocessing as mp
import networkx as nx, matplotlib.pyplot as plt

def red(a, b):
    return (a - 1) % b + 1

n=9
with open(f"nonassociative-affine-triples-{n}.txt", "r") as f:
    unclean_triples = f.read()[2:-2].split('], [')
triples = []
for unclean in unclean_triples:
    t = str.maketrans({',': ' ', '[': ' ', ']': ' '})
    tuple = list(map(int, unclean.translate(t).split()))
    if tuple:
        triples += [triple.BDTriple(tuple)]

counter = 0
for trip in triples:
    counter += 1
    tup = trip.tuple
    edges = []
    G = nx.DiGraph()
    for i in range(n):
        G.add_edge(f"{i+1}", f"{red(i+2, n)}", color='w', weight=0)
    for i in range(n):
        if tup[i] != 0:
            G.add_edge(f"{i+1}", f"{tup[i]}", color='b', weight=2)
    pos = nx.circular_layout(G)
    colors = nx.get_edge_attributes(G,'color').values()
    weights = nx.get_edge_attributes(G,'weight').values()
    nx.draw(G, pos=pos, with_labels=True, edge_color=colors, width=list(weights))
    plt.savefig(f'graph{n}/{counter}-graph.png')
    plt.close()