import mat2
# this is a circular import. to be resolved
from numpy import zeros
import sympy as sp
from copy import deepcopy

def hash(dim, i, j):
    return int(i * dim + j)

def np_to_mat1(dim, np_matrix):
    coef = sp.MutableDenseNDimArray(np_matrix.astype(int))
    return MatrixTensor1(dim, coef)

def to_sparray(dim, coef):
    temp = sp.MutableDenseNDimArray(zeros((dim,)*2).astype(int))
    for i, j in [(x, y) for x in range(dim) for y in range(dim)]:
        temp[i, j] = coef[hash(dim, i, j)]
    return temp

def zero(dim):
    coef = sp.MutableDenseNDimArray(zeros((dim,)*2).astype(int))
    return MatrixTensor1(dim, coef)

def e(dim, i, j):
    coef = sp.MutableDenseNDimArray(zeros((dim,)*2).astype(int))
    coef[i, j] = 1
    return MatrixTensor1(dim, coef)

def gram_schmidt(basis):
        ortho = [basis[0]]
        for i in range(1, len(basis)):
            temp = basis[i]
            for j in range(i):
                temp -= basis[i].pr(ortho[j])
            ortho.append(temp)
        return ortho

class MatrixTensor1:
    def __init__(self, dim, coef):
        if len(coef) == pow(dim, 2):
            self.dim = dim
            self.coef = coef
        else:
            raise Exception("Dimension error")
    
    def __eq__(self, other):
        return self.coef == other.coef

    def __repr__(self):
        dim = self.dim
        coef = self.coef
        sub = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
        result = ""
        for i, j in [(x, y) for x in range(dim) for y in range(dim)]:
            temp = coef[i, j]
            term = f"e{i+1},{j+1}".translate(sub)
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
    
    # Only use this for matrix of GGS type!
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
            return mat2.MatrixTensor2(dim, coef, True)
    
    def pr_to_sln(self):
        dim = self.dim
        coef1 = self.coef
        basis = [e(dim, i + 1, i + 1) - e(dim, i, i) for i in range(dim - 1)]
        ortho_basis = gram_schmidt(basis)
        res_coef = to_sparray(dim, [0] * pow(dim, 2))
        for v in ortho_basis:
            res_coef = res_coef + self.pr(v).coef
        return MatrixTensor1(dim, res_coef)
        coef = sp.MutableDenseNDimArray(zeros((dim,)*2).astype(int))
        for i, j in [(x, y) for x in range(dim) for y in range(dim)]:
            if i != j:
                coef[i, j] += coef1[i, j]
            else:
                if i < dim-1:
                    coef[i+1, i+1] -= sp.Rational(1, 2) * coef1[i, i]
                coef[i, i] += (int(i < dim-1) * sp.Rational(1, 2) 
                               + int(i > 0) * sp.Rational(1, 2)) * coef1[i, i]
                if i > 0:
                    coef[i-1, i-1] -= sp.Rational(1, 2) * coef1[i, i]
        return MatrixTensor1(dim, coef)

    def trace(self):
        tr = 0
        for i in range(self.dim):
            tr += self.coef[i, i]
        return tr
    
    def pr(self, w):
        return sp.Rational((self * w).trace(), (w * w).trace()) * w