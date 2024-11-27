from numpy import zeros
import sympy as sp
from copy import deepcopy

def hash(dim, i, j, k, l, p, q):
        return i * pow(dim, 5) + j * pow(dim, 4) + k * pow(dim, 3) + l * dim * dim + p * dim + q

def hash_human(dim, i, j, k, l, p, q):
    return hash(dim, i-1, j-1, k-1, l-1, p-1, q-1)

def to_sparray(dim, coef):
    temp = sp.MutableDenseNDimArray(zeros((dim,)*6))
    for i, j in [(x, y) for x in range(dim) for y in range(dim)]:
            for k, l in [(x, y) for x in range(dim) for y in range(dim)]:
                for p, q in [(x, y) for x in range(dim) for y in range(dim)]:
                    temp[i, j, k, l, p, q] = coef[hash(dim, i, j, k, l, p, q)]
    return temp

class MatrixTensor3:
    def __init__(self, dim, coef):
        if len(coef) == pow(dim, 6): 
            self.dim = dim
            self.coef = coef
        else:
            raise(ValueError, "Dimension error")

    def __str__(self):
        dim = self.dim
        coef = self.coef
        result = ""
        for i, j in [(x, y) for x in range(dim) for y in range(dim)]:
            for k, l in [(x, y) for x in range(dim) for y in range(dim)]:
                for p, q in [(x, y) for x in range(dim) for y in range(dim)]:
                    temp = coef[i, j, k, l, p, q]
                    sub = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
                    term = f"e{i+1}{j+1}⊗e{k+1}{l+1}⊗e{p+1}{q+1}".translate(sub)
                    if temp != 0:
                        result += "(" + str(temp) + ") * " + term + " + "
        if result == "":
            return "0"
        else:
            return(result[:-3])

    def __add__(self, other):
        if self.dim == other.dim:
            dim = self.dim
            coef1 = self.coef
            coef2 = other.coef
            return MatrixTensor3(dim, coef1 + coef2)
        else: 
            raise(ValueError, "Cannot add matrices of different dimensions")
    
    def __mul__(self, other):
        if self.dim == other.dim:
            dim = self.dim
            coef1 = self.coef
            coef2 = other.coef
            coef = sp.MutableDenseNDimArray(zeros((dim,)*6))
            for i, j in [(x, y) for x in range(dim) for y in range(dim)]:
                for k, l in [(x, y) for x in range(dim) for y in range(dim)]:
                    for p, q in [(x, y) for x in range(dim) for y in range(dim)]:
                        for x, y, z in [(x, y, z) for x in range(dim) for y in range(dim) for z in range(dim)]:
                            coef[i, j, k, l, p, q] += coef1[i, x, k, y, p, z] * coef2[x, j, y, l, z, q]
            return MatrixTensor3(dim, coef)
        else: 
            raise(ValueError, "Cannot multiply matrices of different dimensions")
    
    def __rmul__(self, other):
        return MatrixTensor3(self.dim, other * self.coef)

    def __sub__(self, other):
        return self + (-1) * other

    def comm(self, other):
        return self * other - other * self

    def subs(self, x, y):
        coef1 = self.coef
        dim = self.dim
        coef = sp.MutableDenseNDimArray(zeros((dim,)*6))
        for i, j in [(x, y) for x in range(dim) for y in range(dim)]:
            for k, l in [(x, y) for x in range(dim) for y in range(dim)]:
                for p, q in [(x, y) for x in range(dim) for y in range(dim)]:
                    temp = coef1[i, j, k, l, p, q]
                    if isinstance(temp, sp.Expr):
                        coef[i, j, k, l, p, q] = temp.subs(x, y)
                    else:
                        coef[i, j, k, l, p, q] = temp
        return MatrixTensor3(dim, coef)

    def simplify(self):
        coef = deepcopy(self.coef)
        dim = self.dim
        for i, j in [(x, y) for x in range(dim) for y in range(dim)]:
            for k, l in [(x, y) for x in range(dim) for y in range(dim)]:
                for p, q in [(x, y) for x in range(dim) for y in range(dim)]:
                    temp = coef[i, j, k, l, p, q]
                    if isinstance(temp, sp.Expr):
                        coef[i, j, k, l, p, q] = temp.simplify()
        return MatrixTensor3(dim, coef)