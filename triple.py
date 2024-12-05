from operator import itemgetter
from itertools import groupby
from numpy import zeros
from copy import deepcopy
import sympy as sp
import mat2, mat

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
    N = n + 1
    aux = []
    for i in basis:
        temp = sp.MutableDenseNDimArray(zeros((N,)*2).astype(int))
        temp[(i-1) % N, (i-1) % N] = -sp.Integer(1)
        temp[i % N, i % N] = sp.Integer(1)
        aux.append(mat.MatrixTensor1(N, temp))
    
    res = []
    for i, j in [(x, y) for x in range(n) for y in range(n)]:
        temp = sp.MutableDenseNDimArray(zeros((N,)*2).astype(int))
        temp_mat = mat.MatrixTensor1(N, temp)
        if i <= j:
            temp_mat += sp.Integer(i + 1) * sp.Rational(n - j, n + 1) * aux[j]
        else:
            temp_mat += sp.Integer(j + 1) * sp.Rational(n - i, n + 1) * aux[i]
        print(type(temp_mat.coef[0, 0]))
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
            for i in range(len(basis)):
                if basis[i] not in g1:
                    dual_basis_modified.append(dual_basis_temp[i])
                else:
                    dual_basis_modified.append(sp.Rational(1, 2) * \
                                               (dual_basis_temp[i] + dual_basis_temp[self.T(i)]))

            def aux_inner(triple, i, j):
                k = triple.T(i)
                l = triple.T(j)
                return 2 * (int(i == j) + int(k == j)) \
                    - 2 * int(l != 0) * (int(i == l) + int(k == l)) \
                    - int(i+1 == j) - int(i == j+1) - int(k+1 == j) - int(k == j+1) \
                    + int(l != 0) * (int(i+1 == l) + int(i == l+1) + int(k+1 == l) + int(k == l+1))

            
            res = mat2.MatrixTensor2(n, sp.MutableDenseNDimArray(zeros((n,)*4).astype(int)))
            for i, j in [(x, y) for x in range(len(basis)) for y in range(len(basis))]:
                if basis[i] in g1:
                    res += sp.Rational(1, 2) * aux_inner(self, basis[i], basis[j]) * \
                        dual_basis_modified[i].tensor(dual_basis_modified[j])
                elif basis[j] in g1:
                    res -= sp.Rational(1, 2) * aux_inner(self, basis[j], basis[i]) * \
                        dual_basis_modified[i].tensor(dual_basis_modified[j])
            return res + sp.Rational(1, 2) * mat2.casimir(n)
    
    # Here x is a Sympy symbol, and r0 is the continuous datum we choose.
    # Not including the computation of r0 in this method prevents us from repeatedly calculating r0.
    def to_trigonometric_solution(self, r0, x) -> mat2.MatrixTensor2:
        n = self.n
        coef = mat2.to_sparray([0] * pow(n, 4))
        for m in range(1, n+1):
            for i, j in [(x, y) for x in range(n) for y in range(n)]:
                if (i - j) % n == m:
                    coef[i, j, j, i] += sp.exp(m * x / n) / (sp.exp(x) - 1)
        return mat2.MatrixTensor2(coef) + r0