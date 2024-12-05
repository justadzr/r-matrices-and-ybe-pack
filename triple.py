from operator import itemgetter
from itertools import groupby
from numpy import zeros
from copy import deepcopy
import sympy as sp
import mat2, mat1
import numpy as np

# Write a BD triple of the affine untwisted sl(n) as a list of length n, indicating the image
# of \alpha_i under the transformation T: for instance
# [a_1, ..., a_n]
# means T(\alpha_1)=\alpha_{a_n}. Note if a_i is then \alpha_i is not in \Gamma_1.
# From now let 1, ..., n represents the roots of the affine untwisted sl(n). 
# Write 0 if the ith root is not in the domain of the transformation.

# The code here is not optimal at all: nested loops and unnecessary if/else are everywhere.
# But usually we don't need large triples (at least our matrix computation can't be too hard)
# since otherwise Sympy.simplify() would be too slow. So I just leave it as it is.

def left_end(n, connected):
    for i in range(len(connected)):
        if (connected[i] - 2 + n) % n + 1 not in connected:
            return connected[i]

def right_end(n, connected):
    for i in range(len(connected)):
        if connected[i] % n + 1 not in connected:
            return connected[i]

# The input basis here is a list of n=N-1 distinct integers from 1 to N representing
# a basis of the Cartan algebra of sl(N)
def dual_basis(basis):
    n = len(basis)
    # print(basis)
    N = n + 1
    aux = []
    for i in basis:
        temp = sp.MutableDenseNDimArray(zeros((N,)*2).astype(int))
        temp[(i-1) % N, (i-1) % N] = -sp.Integer(1)
        temp[i % N, i % N] = sp.Integer(1)
        aux.append(mat1.MatrixTensor1(N, temp))
    
    system = sp.zeros(n, n)
    for i, j in [(x, y) for x in range(n) for y in range(n)]:
        system[i, j] = 2 * int(i == j) - int(basis[j] == basis[i] % N + 1) \
            - int(basis[j] == (basis[i] -2 + N) % N + 1)
    inverse = system.inv()
    res = []
    for i in range(n):
        temp = sp.MutableDenseNDimArray(zeros((N,)*2).astype(int))
        temp_mat = mat1.MatrixTensor1(N, temp)
        for j in range(n):
            temp_mat += inverse[i, j] * aux[j]
        res.append(temp_mat)
    return res

class BDTriple:
    def __init__(self, triple):
        self.n = len(triple)
        self.g1 = [(x+1) for x in range(self.n) if triple[x] != 0]
        self.g2 = [triple[x-1] for x in self.g1]

    def valid(self) -> bool:
        if len(set(self.g2)) != len(self.g1):
            print("Not bijective")
            return False

        n = self.n
        for i in range(len(self.g1)):
            if self.g1[i] % n + 1 in self.g1:
                ind = self.g1.index(self.g1[i] % n + 1)
                temp = abs(self.g2[i] - self.g2[ind])
                if temp != 1 and temp != n - 1:
                    print("Not orthogonal")
                    return False
            if (self.g1[i] - 2 + n) % n + 1 in self.g1:
                ind = self.g1.index((self.g1[i] - 2 + n) % n + 1)
                temp = abs(self.g2[i] - self.g2[ind])
                if temp != 1 and temp != n - 1:
                    print("Not orthogonal")
                    return False
            if self.g2[i] % n + 1 in self.g2:
                ind = self.g2.index(self.g2[i] % n + 1)
                temp = abs(self.g1[i] - self.g1[ind])
                if temp != 1 and temp != n - 1:
                    print("Not orthogonal")
                    return False
            if (self.g2[i] - 2 + n) % n + 1 in self.g2:
                ind = self.g2.index((self.g2[i] - 2 + n) % n + 1)
                temp = abs(self.g1[i] - self.g1[ind])
                if temp != 1 and temp != n - 1:
                    print("Not orthogonal")
                    return False

        for i in range(len(self.g1)):
            temp = self.g1[i]
            for k in range(self.n + 1):
                temp = self.T(temp)
                if temp == 0:
                    break
            if k >= self.n:
                print("Not nilpotent")
                return False
        return True
    
    def T(self, i):
        if i not in self.g1:
            return 0
        else:
            return self.g2[self.g1.index(i)]

    def __repr__(self):
        res = ""
        sub = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
        for i in range(len(self.g1)):
            res += f"\u03B1{self.g1[i]} -> \u03B1{self.g2[i]}\n".translate(sub)
        return res

    def connected_components_aux(self, G):
        n = self.n
        res = [None] * n
        m = len(G)
        key_1, key_n = None, None
        for k, g in groupby(enumerate(G), lambda x : (x[0] - x[1]) % m):
            temp = list(map(itemgetter(1), g))
            if 1 in temp:
                key_1 = k
            if n in temp:
                key_n = k
            if res[k] is None:
                res[k] = temp
            else:
                res[k] += temp
        if key_1 is not None and key_n is not None and key_1 != key_n:
            res[key_1] += res[key_n]
            res[key_n] = None
        return [x for x in res if x is not None]
    
    def connected_components(self):
        return self.connected_components_aux(self.g1)
    
    def connected_components_img(self):
        components = self.connected_components()
        unclean = self.connected_components_aux(sorted(self.g2))
        temp = []
        if len(components) != len(unclean):
            print(unclean)
            raise Exception("Triple not valid: not orthogonal.")
        else: 
            for i in components:
                for j in unclean:
                    if self.T(i[0]) in j:
                        temp.append(j)
                        break
            return temp

    def associative(self) -> bool:
        if not self.valid():
            raise Exception("Triple not valid")
        else:
            components = self.connected_components()
            components_img = self.connected_components_img()
            n = self.n
            for i in range(len(components)):
                length = len(components[i])
                if 1 < length < n: 
                    # It is impossible for the length to be n. But I will leave it be.
                    l1 = left_end(n, components[i])
                    l2 = left_end(n, components_img[i])
                    if l2 != self.T(l1):
                        print(f"The orientation of" + \
                              " the connected component {components[i]} is reversed")
                        return False
            return True
    
    # We first choose s := r0 - Casimir/2 via Schedler's method and add a half of the Casimir.
    # Please only use this method once for each triple!
    def choose_r0(self) -> mat2.MatrixTensor2:
        if not self.valid():
            raise Exception("Triple not valid")
        else:
            n = self.n
            g1 = self.g1
            g2 = self.g2
            excl = -1
            for i in range(1, n+1):
                if i not in g1:
                    excl = i
                    break
            basis = [x for x in range(1, n+1) if x != excl]
            dual_basis_temp = dual_basis(basis)
            dual_basis_modified = []

            T_mat = sp.zeros(len(basis), len(basis))
            for i in range(len(basis)):
                if basis[i] in g1:
                    j = self.T(basis[i])
                    if j not in basis:
                        for p in range(len(basis)):
                            T_mat[p, i] = -1
                    else:
                        T_mat[basis.index(j), i] = 1
            transform = (sp.eye(len(basis)) - T_mat).inv()

            for i in range(len(basis)):
                temp = mat1.MatrixTensor1(n, sp.MutableDenseNDimArray(zeros((n,)*2).astype(int)))
                for j in range(len(basis)):
                    temp += transform[i, j] * dual_basis_temp[j]
                dual_basis_modified.append(temp)

            def aux_inner(triple, i, j):
                k = triple.T(i)
                l = triple.T(j)
                return 2 * (int(i == j) + int(k == j)) \
                    - 2 * int(l != 0) * (int(i == l) + int(k == l)) \
                    - int(i%n+1 == j) - int(i == j%n+1) - int(k%n+1 == j) - int(k == j%n+1) \
                    + int(l != 0) * (int(i%n+1 == l) + int(i == l%n+1) + \
                                     int(k%n+1 == l) + int(k == l%n+1))
            
            s = mat2.MatrixTensor2(n, sp.MutableDenseNDimArray(zeros((n,)*4).astype(int)))
            for i, j in [(x, y) for x in range(len(basis)) for y in range(len(basis))]:
                if basis[i] in g1:
                    s += sp.Rational(1, 2) * aux_inner(self, basis[i], basis[j]) * \
                        dual_basis_modified[i].tensor(dual_basis_modified[j])
                elif basis[j] in g1:
                    s -= sp.Rational(1, 2) * aux_inner(self, basis[j], basis[i]) * \
                        dual_basis_modified[i].tensor(dual_basis_modified[j])
            return s + sp.Rational(1, 2) * mat2.casimir(n)