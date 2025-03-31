import sympy as sp, triple, ybe, multiprocessing as mp
import matplotlib

n=10
with open(f"nonassociative-affine-triples-{n}.txt", "r") as f:
    unclean_triples = f.read()[2:-2].split('], [')
triples = []
for unclean in unclean_triples:
    t = str.maketrans({',': ' ', '[': ' ', ']': ' '})
    tuple = list(map(int, unclean.translate(t).split()))
    if tuple:
        triples += [triple.BDTriple(tuple)]