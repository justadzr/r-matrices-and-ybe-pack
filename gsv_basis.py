from triple import BDTriple
from mat2 import MatrixTensor2, zero
from mat1 import np_to_mat1, MatrixTensor1
import numpy as np
import jax

def red(a, b):
    return (a - 1) % b + 1

def standard_form_mat2(m1 : MatrixTensor2, m2 : MatrixTensor2):
    n = m1.dim
    coef1 = m1.coef
    coef2 = m2.coef
    res = 0
    for (i, j) in [(x, y) for x in range(n) for y in range(n)]:
        for (k, l) in [(x, y) for x in range(n) for y in range(n)]:
            res = coef1[i, j, k, l] * coef2[j, i, l, k] + res
    return res

def DR_standard_form(f):
    return 0

def DL_standard_form():
    return 0

def gsv_omega(dim, r1, r2, fi, fj, bracket):
    # search for a matrix g such that f_i and f_j are both nonzero on g
    def find_nonzero(n, fi, fj, trials=1000, tol=1e-12):
        for _ in range(trials):
            g_np = np.random.randn(n, n)
            det = np.linalg.det(g_np)
            if abs(det) > tol:
                g = np_to_mat1(n, g_np / det ** (1/n))
                vali, valj = fi(g), fj(g)
                if abs(vali) > tol and abs(valj) > tol:  # avoid numerical noise near 0
                    return g
    g = find_nonzero(dim, fi, fj)
    return bracket(dim, r1, r2, fi, fj, g) / (fi(g) * fj(g))

# def frozen(i, j):
#     return i == 0 or j == 0

# def gln_frozen(i, j):
#     return i == 0 and j == 0

def runs(n, g, components):
    runs = []
    dual_runs = []
    components_ordered = sorted(components, key=lambda x: x[0])
    for comp in components_ordered:
        temp = comp + [red(comp[-1] + 1, n)]
        runs.append(temp)

    for i in range(1, n+1):
        if i not in g:
            print(i, runs)
            if len(runs) == 1:
                if i < runs[0][0]:
                    runs = [[i]] + runs
                if i > runs[0][-1]:
                    runs = runs + [[i]]
            for ind in range(len(runs)):
                if runs[ind][-1] < i:
                    if ind + 1 < len(runs) and i < runs[ind+1][0]:
                        runs = runs[:ind+1] + [[i]] + runs[ind+1:]
                    elif ind + 1 == len(runs) and runs[ind][-1] < i:
                        runs = runs + [[i]]
    for run in runs:
        dual_runs.append([x for x in range(n - run[-1] + 1, n - run[0] + 2)])
    return runs, dual_runs

def X_runs(triple: BDTriple):
    return runs(triple.n, triple.g1, triple.connected_components())

def Y_runs(triple: BDTriple):
    return runs(triple.n, triple.g2, triple.connected_components_img())

def Y_runs_corrected(triple: BDTriple):
    x_runs = X_runs(triple)
    y_runs = Y_runs(triple)

def gsv_quiver(triple: BDTriple) -> MatrixTensor2:
    dim = triple.n
    coef = zero(dim).coef
    for i in range(dim):
        for j in range(dim):
            if i < dim - 1 and j < dim - 1:
                coef[i, j, i+1, j+1] += 1
                coef[i+1, j+1, i, j] -= 1
            if i != 0 and j != 0:
                coef[i, j, i-1, j] += 1
                coef[i-1, j, i, j] -= 1
                coef[i, j, i, j-1] += 1
                coef[i, j-1, i, j] -= 1
    return MatrixTensor2(dim, coef)

def sklyanin(dim, r1, r2, f_i, f_j, g):
    return 0

class GSV(BDTriple):
    def __init__(self, tuple, **kwarg):
        super().__init__(tuple, **kwarg)
