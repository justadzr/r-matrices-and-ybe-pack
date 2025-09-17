import sympy as sp
from copy import deepcopy
from numpy import zeros
from time import time
from mat1 import e

def hash(dim, i, j, k, l, p, q):
        return i * pow(dim, 5) + j * pow(dim, 4) + k * pow(dim, 3) + l * dim * dim + p * dim + q

def hash_human(dim, i, j, k, l, p, q):
    return hash(dim, i-1, j-1, k-1, l-1, p-1, q-1)

def to_sparray(dim, coef):
    temp = sp.MutableDenseNDimArray(zeros((dim,)*6).astype(int))
    for i, j in [(x, y) for x in range(dim) for y in range(dim)]:
            for k, l in [(x, y) for x in range(dim) for y in range(dim)]:
                for p, q in [(x, y) for x in range(dim) for y in range(dim)]:
                    temp[i, j, k, l, p, q] = coef[hash(dim, i, j, k, l, p, q)]
    return temp

def zero(dim):
    coef = sp.MutableDenseNDimArray(zeros((dim,)*6).astype(int))
    return MatrixTensor3(dim, coef, True)

class MatrixTensor3:
    def __init__(self, dim: int, coef: sp.MutableDenseNDimArray, of_ggs_type: bool):
        if len(coef) == pow(dim, 6): 
            self.dim = dim
            self.coef = coef
            self.ggs = of_ggs_type
        else:
            raise(ValueError, "Dimension error")

    def __repr__(self):
        dim = self.dim
        coef = self.coef
        result = ""
        for i, j in [(x, y) for x in range(dim) for y in range(dim)]:
            for k, l in [(x, y) for x in range(dim) for y in range(dim)]:
                for p, q in [(x, y) for x in range(dim) for y in range(dim)]:
                    temp = coef[i, j, k, l, p, q]
                    sub = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
                    term = f"e{i+1},{j+1}⊗e{k+1},{l+1}⊗e{p+1},{q+1}".translate(sub)
                    if temp != 0:
                        result += "(" + str(sp.nsimplify(temp)) + ") * " + term + " + "
        if result == "":
            return "0"
        else:
            return(result[:-3])

    def __eq__(self, other):
        return self.coef == other.coef #and self.ggs == other.ggs

    def __add__(self, other):
        if self.dim == other.dim:
            dim = self.dim
            coef1 = self.coef
            coef2 = other.coef
            return MatrixTensor3(dim, coef1 + coef2, self.ggs and other.ggs)
        else: 
            raise(ValueError, "Cannot add matrices of different dimensions")
    
    # Please do not multiply MatrixTensor3
    def __mul__(self, other):
        start = time()
        if self.dim == other.dim:
            dim = self.dim
            coef1 = self.coef
            coef2 = other.coef
            coef = sp.MutableDenseNDimArray(zeros((dim,)*6).astype(int))
            # if self.ggs and other.ggs:
            #     counter = 0
            #     for i, j in [(x, y) for x in range(dim) for y in range(dim)]:
            #         for k, l in [(x, y) for x in range(dim) for y in range(dim)]:
            #             for p, q in [(x, y) for x in range(dim) for y in range(dim)]:
            #                 for x, y, z in [(x, y, z) for x in range(dim) for y in range(dim) for z in range(dim)]:
            #                     coef[i, j, k, l, p, q] += coef1[i, x, k, y, p, z] * coef2[x, j, y, l, z, q]
            #                     counter += 1
            # else:
            counter = 0
            for i, j in [(x, y) for x in range(dim) for y in range(dim)]:
                for k, l in [(x, y) for x in range(dim) for y in range(dim)]:
                    for p, q in [(x, y) for x in range(dim) for y in range(dim)]:
                        for x, y, z in [(x, y, z) for x in range(dim) for y in range(dim) for z in range(dim)]:
                            coef[i, j, k, l, p, q] += coef1[i, x, k, y, p, z] * coef2[x, j, y, l, z, q]
                            counter += 1
            print(f"It took us {time() - start} seconds to do {counter} multiplications.")
            return MatrixTensor3(dim, coef, self.ggs and other.ggs)
        else: 
            raise(ValueError, "Cannot multiply matrices of different dimensions")
    
    def __rmul__(self, other):
        return MatrixTensor3(self.dim, other * self.coef, self.ggs)
    
    def __sub__(self, other):
        return self + (-1) * other

    def __add__(self, other):
        if self.dim == other.dim:
            dim = self.dim
            coef1 = self.coef
            coef2 = other.coef
            return MatrixTensor3(dim, coef1 + coef2, self.ggs and other.ggs)
        else: 
            raise(ValueError, "Cannot add matrices of different dimensions")
    
    def __neg__(self):
        return (-1) * self

    def comm(self, other):
        return self * other - other * self

    def subs(self, x, y):
        coef1 = self.coef
        dim = self.dim
        coef = sp.MutableDenseNDimArray(zeros((dim,)*6).astype(int))
        for i, j in [(x, y) for x in range(dim) for y in range(dim)]:
            for k, l in [(x, y) for x in range(dim) for y in range(dim)]:
                for p, q in [(x, y) for x in range(dim) for y in range(dim)]:
                    temp = coef1[i, j, k, l, p, q]
                    if isinstance(temp, sp.Expr):
                        coef[i, j, k, l, p, q] = temp.subs(x, y)
                    else:
                        coef[i, j, k, l, p, q] = temp
        return MatrixTensor3(dim, coef, self.ggs)

    def simplify(self):
        start = time()
        coef = deepcopy(self.coef)
        dim = self.dim
        for i, j in [(x, y) for x in range(dim) for y in range(dim)]:
            for k, l in [(x, y) for x in range(dim) for y in range(dim)]:
                for p, q in [(x, y) for x in range(dim) for y in range(dim)]:
                    temp = coef[i, j, k, l, p, q]
                    if isinstance(temp, sp.Expr):
                        coef[i, j, k, l, p, q] = temp.simplify()
        return MatrixTensor3(dim, coef, self.ggs)
    
    def simplify_rat(self):
        start = time()
        coef = deepcopy(self.coef)
        dim = self.dim
        for i, j in [(x, y) for x in range(dim) for y in range(dim)]:
            for k, l in [(x, y) for x in range(dim) for y in range(dim)]:
                for p, q in [(x, y) for x in range(dim) for y in range(dim)]:
                    temp = coef[i, j, k, l, p, q]
                    if isinstance(temp, sp.Expr):
                        coef[i, j, k, l, p, q] = temp.ratsimp()
        return MatrixTensor3(dim, coef, self.ggs)
    
    def check_of_ggs_type(self):
        dim = self.dim
        coef = self.coef
        for i, j in [(x, y) for x in range(dim) for y in range(dim)]:
            for k, l in [(x, y) for x in range(dim) for y in range(dim)]:
                for p, q in [(x, y) for x in range(dim) for y in range(dim)]:
                    if q != (i + k + p - j - l) % dim and coef[i, j, k, l, p, q] != 0:
                        return MatrixTensor3(dim, coef, False)
        return MatrixTensor3(dim, coef, True)

    def pr_to_sln(self):
        dim = self.dim
        coef = self.coef
        res = zero(dim)
        for i, j in [(x, y) for x in range(dim) for y in range(dim)]:
            for k, l in [(x, y) for x in range(dim) for y in range(dim)]:
                for p, q in [(x, y) for x in range(dim) for y in range(dim)]:
                    res += coef[i, j, k, l, p, q] * \
                        (e(dim, i, j).pr_to_sln().tensor(e(dim, k, l).pr_to_sln())).tensor(e(dim, p, q).pr_to_sln())

        return res
    
    def subs(self, x, y):
        dim = self.dim
        coef = self.coef
        coef1 = deepcopy(coef)
        for i, j in [(x, y) for x in range(dim) for y in range(dim)]:
            for k, l in [(x, y) for x in range(dim) for y in range(dim)]:
                for p, q in [(x, y) for x in range(dim) for y in range(dim)]:
                    coef1[i, j, k, l, p, q].subs(x, y)
        return MatrixTensor3(dim, coef1, self.ggs)