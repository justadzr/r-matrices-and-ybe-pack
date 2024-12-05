import mat2
from numpy import zeros
import sympy as sp
from copy import deepcopy

def hash(dim, i, j):
    return int(i * dim + j)

def to_sparray(dim, coef):
    temp = sp.MutableDenseNDimArray(zeros((dim,)*4).astype(int))
    for i, j in [(x, y) for x in range(dim) for y in range(dim)]:
        temp[i, j] = coef[hash(dim, i, j)]
    return temp

class MatrixTensor1:
    def __init__(self, dim, coef):
        if len(coef) == pow(dim, 2):
            self.dim = dim
            self.coef = coef
        else:
            raise Exception("Dimension error")

    def __repr__(self):
        dim = self.dim
        coef = self.coef
        sub = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
        result = ""
        for i, j in [(x, y) for x in range(dim) for y in range(dim)]:
            temp = coef[i, j]
            term = f"e{i+1}{j+1}".translate(sub)
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
            return MatrixTensor1(dim, coef1 + coef2)
        else: 
            raise(ValueError, "Cannot add matrices of different dimensions")

    def __rmul__(self, other):
        return MatrixTensor1(self.dim, other * self.coef)

    def __sub__(self, other):
        return self + (-1) * other
    
    def __mul__(self, other):
        if self.dim == other.dim:
            dim = self.dim
            coef1 = self.coef
            coef2 = other.coef
            coef = sp.MutableDenseNDimArray(zeros((dim,)*2).astype(int))
            for i, j in [(x, y) for x in range(dim) for y in range(dim)]:
                for p, q in [(x, y) for x in range(dim) for y in range(dim)]:
                    coef[i, j] += coef1[i, p] * coef2[q, j]
            return MatrixTensor1(dim, coef)
        else: 
            raise(ValueError, "Cannot multiply matrices of different dimensions")
    
    def comm(self, other):
        return self * other - other * self

    def subs(self, x, y):
        coef1 = self.coef
        dim = self.dim
        coef = sp.MutableDenseNDimArray(zeros((dim,)*2).astype(int))
        for i, j in [(x, y) for x in range(dim) for y in range(dim)]:
            temp = coef1[i, j]
            if isinstance(temp, sp.Expr):
                coef[i, j] = temp.subs(x, y)
            else:
                coef[i, j] = temp
        return MatrixTensor1(dim, coef)
    
    def simplify(self):
        coef = deepcopy(self.coef)
        dim = self.dim
        for i, j in [(x, y) for x in range(dim) for y in range(dim)]:
            temp = coef[i, j]
            if isinstance(temp, sp.Expr):
                coef[i, j] = temp.simplify()
        return MatrixTensor1(dim, coef)
    
    def tensor(self, other):
        if self.dim != other.dim:
            raise(ValueError, "Cannot tensor matrices of different dimensions")
        else: 
            coef1 = self.coef
            coef2 = other.coef
            dim = self.dim
            coef = sp.MutableDenseNDimArray(zeros((dim,)*4).astype(int))
            for i, j in [(x, y) for x in range(dim) for y in range(dim)]:
                for k, l in [(x, y) for x in range(dim) for y in range(dim)]:
                    coef[i, j, k, l] = coef1[i, j] * coef2[k, l]
            return mat2.MatrixTensor2(dim, coef)