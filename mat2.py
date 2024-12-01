import mat3
from numpy import zeros
import sympy as sp
from copy import deepcopy

def hash(dim, i, j, k, l):
    return int(i * pow(dim, 3) + j * dim * dim + k * dim + l)

def hash_human(dim, i, j, k, l):
    return hash(dim, i-1, j-1, k-1, l-1)

def to_sparray(dim, coef):
    temp = sp.MutableDenseNDimArray(zeros((dim,)*4))
    for i, j in [(x, y) for x in range(dim) for y in range(dim)]:
            for k, l in [(x, y) for x in range(dim) for y in range(dim)]:
                temp[i, j, k, l] = coef[hash(dim, i, j, k, l)]
    return temp

def casimir(dim):
        coef = sp.MutableDenseNDimArray(zeros((dim,)*4))
        for i, j in [(x, y) for x in range(dim) for y in range(dim)]:
            coef[i, j, j, i] = 1
        return MatrixTensor2(dim, coef)

class MatrixTensor2:
    def __init__(self, dim, coef):
        if len(coef) == pow(dim, 4): 
            self.dim = dim
            self.coef = coef
        else:
            raise(ValueError, "Dimension error")

    def __str__(self):
        dim = self.dim
        coef = self.coef
        sub = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
        result = ""
        for i, j in [(x, y) for x in range(dim) for y in range(dim)]:
            for k, l in [(x, y) for x in range(dim) for y in range(dim)]:
                temp = coef[i, j, k, l]
                term = f"e{i+1}{j+1}⊗e{k+1}{l+1}".translate(sub)
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
            return MatrixTensor2(dim, coef1 + coef2)
        else: 
            raise(ValueError, "Cannot add matrices of different dimensions")

    def __rmul__(self, other):
        return MatrixTensor2(self.dim, other * self.coef)

    def __sub__(self, other):
        return self + (-1) * other
    
    def __mul__(self, other):
        if self.dim == other.dim:
            dim = self.dim
            coef1 = self.coef
            coef2 = other.coef
            coef = sp.MutableDenseNDimArray(zeros((dim,)*4))
            for i, j in [(x, y) for x in range(dim) for y in range(dim)]:
                for k, l in [(x, y) for x in range(dim) for y in range(dim)]:
                    for p, q in [(x, y) for x in range(dim) for y in range(dim)]:
                        coef[i, j, k, l] += coef1[i, p, k, q] * coef2[p, j, q, l]
            return MatrixTensor2(dim, coef)
        else: 
            raise(ValueError, "Cannot multiply matrices of different dimensions")
    
    def swap(self):
        coef1 = self.coef
        dim = self.dim
        coef = sp.MutableDenseNDimArray(zeros((dim,)*4))
        for i, j in [(x, y) for x in range(dim) for y in range(dim)]:
                for k, l in [(x, y) for x in range(dim) for y in range(dim)]:
                    coef[i, j, k, l] = coef1[k, l, i, j]
        return MatrixTensor2(dim, coef)
    
    def flip(self):
        return casimir(self.dim) * self
    
    def to_matrixtensor3_12(self):
        coef1 = self.coef
        dim = self.dim
        coef = sp.MutableDenseNDimArray(zeros((dim,)*6))
        for i, j in [(x, y) for x in range(dim) for y in range(dim)]:
                for k, l in [(x, y) for x in range(dim) for y in range(dim)]:
                    for p in range(dim):
                        coef[i, j, k, l, p, p] = coef1[i, j, k, l]
        return mat3.MatrixTensor3(dim, coef)

    def to_matrixtensor3_23(self):
        coef1 = self.coef
        dim = self.dim
        coef = sp.MutableDenseNDimArray(zeros((dim,)*6))
        for i, j in [(x, y) for x in range(dim) for y in range(dim)]:
                for k, l in [(x, y) for x in range(dim) for y in range(dim)]:
                    for p in range(dim):
                        coef[p, p, i, j, k, l] = coef1[i, j, k, l]
        return mat3.MatrixTensor3(dim, coef)
        
    def to_matrixtensor3_13(self):
        coef1 = self.coef
        dim = self.dim
        coef = sp.MutableDenseNDimArray(zeros((dim,)*6))
        for i, j in [(x, y) for x in range(dim) for y in range(dim)]:
                for k, l in [(x, y) for x in range(dim) for y in range(dim)]:
                    for p in range(dim):
                        coef[i, j, p, p, k, l] = coef1[i, j, k, l]
        return mat3.MatrixTensor3(dim, coef)
    
    def comm(self, other):
        return self * other - other * self

    def subs(self, x, y):
        coef1 = self.coef
        dim = self.dim
        coef = sp.MutableDenseNDimArray(zeros((dim,)*4))
        for i, j in [(x, y) for x in range(dim) for y in range(dim)]:
            for k, l in [(x, y) for x in range(dim) for y in range(dim)]:
                temp = coef1[i, j, k, l]
                if isinstance(temp, sp.Expr):
                    coef[i, j, k, l] = temp.subs(x, y)
                else:
                    coef[i, j, k, l] = temp
        return MatrixTensor2(dim, coef)
    
    def simplify(self):
        coef = deepcopy(self.coef)
        dim = self.dim
        for i, j in [(x, y) for x in range(dim) for y in range(dim)]:
            for k, l in [(x, y) for x in range(dim) for y in range(dim)]:
                temp = coef[i, j, k, l]
                if isinstance(temp, sp.Expr):
                    coef[i, j, k, l] = temp.simplify()
        return MatrixTensor2(dim, coef)
    
    def exp(self):
        return 0