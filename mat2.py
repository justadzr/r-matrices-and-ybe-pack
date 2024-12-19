import mat3, mat1, sympy as sp
from numpy import zeros
from copy import deepcopy

# import time

def hash(dim, i, j, k, l):
    return int(i * pow(dim, 3) + j * dim * dim + k * dim + l)

def hash_human(dim, i, j, k, l):
    return hash(dim, i-1, j-1, k-1, l-1)

def to_sparray(dim, coef):
    temp = sp.MutableDenseNDimArray(zeros((dim,)*4).astype(int))
    for i, j in [(x, y) for x in range(dim) for y in range(dim)]:
            for k, l in [(x, y) for x in range(dim) for y in range(dim)]:
                temp[i, j, k, l] = coef[hash(dim, i, j, k, l)]
    return temp

def casimir(dim):
    coef = sp.MutableDenseNDimArray(zeros((dim,)*4).astype(int))
    res = MatrixTensor2(dim, coef, False)
    for i, j in [(x, y) for x in range(dim) for y in range(dim)]:
        if i < j:
            temp = sp.MutableDenseNDimArray(zeros((dim,)*2).astype(int))
            temp[i, i] = 1
            temp[j, j] = -1
            res += sp.Rational(1, dim) * \
                mat1.MatrixTensor1(dim, temp).tensor(mat1.MatrixTensor1(dim, temp))
    return MatrixTensor2(dim, res.coef, True)

def identity(dim):
    coef = sp.MutableDenseNDimArray(zeros((dim,)*2).astype(int))
    for i in range(dim):
        coef[i, i] = 1
    return mat1.MatrixTensor1(coef).tensor(mat1.MatrixTensor1(coef))

def zero(dim):
    coef = sp.MutableDenseNDimArray(zeros((dim,)*4).astype(int))
    return MatrixTensor2(dim, coef, True)

class MatrixTensor2:
    def __init__(self, dim: int, coef: sp.MutableDenseNDimArray, of_ggs_type: bool):
        if len(coef) == pow(dim, 4): 
            self.dim = dim
            self.coef = coef
            self.ggs = of_ggs_type
        else:
            raise(ValueError, "Dimension error")
    
    def __eq__(self, other):
        return self.coef == other.coef #and self.ggs == other.ggs

    def __repr__(self):
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
            return MatrixTensor2(dim, coef1 + coef2, self.ggs and other.ggs)
        else: 
            raise(ValueError, "Cannot add matrices of different dimensions")

    def __rmul__(self, other):
        return MatrixTensor2(self.dim, other * self.coef, self.ggs)

    def __sub__(self, other):
        return self + (-1) * other
    
    def __mul__(self, other):
        # start_time = time.time()
        if self.dim == other.dim:
            dim = self.dim
            coef1 = self.coef
            coef2 = other.coef
            coef = sp.MutableDenseNDimArray(zeros((dim,)*4).astype(int))
            if self.ggs and other.ggs:
                for i, j in [(x, y) for x in range(dim) for y in range(dim)]:
                    for k in range(dim):
                        l = (i + k - j) % dim
                        for m in range(dim):
                            l_m = (i + k - m) % dim
                            coef[i, j, k, l] += coef1[i, m, k, l_m] * coef2[m, j, l_m, l]
            else:
                for i, j in [(x, y) for x in range(dim) for y in range(dim)]:
                    for k, l in [(x, y) for x in range(dim) for y in range(dim)]:
                        for p, q in [(x, y) for x in range(dim) for y in range(dim)]:
                            coef[i, j, k, l] += coef1[i, p, k, q] * coef2[p, j, q, l]
            # delta_t = time.time() - start_time
            # print(f"A multiplication took: {delta_t} seconds")
            return MatrixTensor2(dim, coef, self.ggs and other.ggs)
        else: 
            raise(ValueError, "Cannot multiply matrices of different dimensions")
    
    def swap(self):
        coef1 = self.coef
        dim = self.dim
        coef = sp.MutableDenseNDimArray(zeros((dim,)*4).astype(int))
        for i, j in [(x, y) for x in range(dim) for y in range(dim)]:
                for k, l in [(x, y) for x in range(dim) for y in range(dim)]:
                    coef[i, j, k, l] = coef1[k, l, i, j]
        return MatrixTensor2(dim, coef, self.ggs)
    
    def flip(self):
        return casimir(self.dim) * self
    
    def to_matrixtensor3_12(self):
        coef1 = self.coef
        dim = self.dim
        coef = sp.MutableDenseNDimArray(zeros((dim,)*6).astype(int))
        for i, j in [(x, y) for x in range(dim) for y in range(dim)]:
                for k, l in [(x, y) for x in range(dim) for y in range(dim)]:
                    for p in range(dim):
                        coef[i, j, k, l, p, p] = coef1[i, j, k, l]
        return mat3.MatrixTensor3(dim, coef, self.ggs)

    def to_matrixtensor3_23(self):
        coef1 = self.coef
        dim = self.dim
        coef = sp.MutableDenseNDimArray(zeros((dim,)*6).astype(int))
        for i, j in [(x, y) for x in range(dim) for y in range(dim)]:
                for k, l in [(x, y) for x in range(dim) for y in range(dim)]:
                    for p in range(dim):
                        coef[p, p, i, j, k, l] = coef1[i, j, k, l]
        return mat3.MatrixTensor3(dim, coef, self.ggs)
        
    def to_matrixtensor3_13(self):
        coef1 = self.coef
        dim = self.dim
        coef = sp.MutableDenseNDimArray(zeros((dim,)*6).astype(int))
        for i, j in [(x, y) for x in range(dim) for y in range(dim)]:
                for k, l in [(x, y) for x in range(dim) for y in range(dim)]:
                    for p in range(dim):
                        coef[i, j, p, p, k, l] = coef1[i, j, k, l]
        return mat3.MatrixTensor3(dim, coef, self.ggs)
    
    def comm(self, other):
        return self * other - other * self

    def subs(self, x, y):
        coef1 = self.coef
        dim = self.dim
        coef = sp.MutableDenseNDimArray(zeros((dim,)*4).astype(int))
        for i, j in [(x, y) for x in range(dim) for y in range(dim)]:
            for k, l in [(x, y) for x in range(dim) for y in range(dim)]:
                temp = coef1[i, j, k, l]
                if isinstance(temp, sp.Expr):
                    coef[i, j, k, l] = temp.subs(x, y)
                else:
                    coef[i, j, k, l] = temp
        return MatrixTensor2(dim, coef, self.ggs)
    
    def identity(self):
        dim = self.dim

    def simplify(self):
        coef = deepcopy(self.coef)
        dim = self.dim
        for i, j in [(x, y) for x in range(dim) for y in range(dim)]:
            for k, l in [(x, y) for x in range(dim) for y in range(dim)]:
                temp = coef[i, j, k, l]
                if isinstance(temp, sp.Expr):
                    coef[i, j, k, l] = temp.simplify()
        return MatrixTensor2(dim, coef, self.ggs)
    
    # The input hbar here is the Sympy.symbol parameter used to define q
    # The input half is a boolean constant indicating if q=e^\hbar or q=e^{\hbar/2}
    # Note the input gets automatically projected to h\otimes h,
    # which is where the exponential makes sense in U_q(g)
    def exp(self, hbar: sp.Symbol, q_exp_half: bool):
        coef1 = self.coef
        dim = self.dim
        res = identity(dim)
        if q_exp_half:
            const = sp.Rational(1, 2)
        else:
            const = 1
        for i, k in [(x, y) for x in range(dim) for y in range(dim)]:
            temp = sp.MutableDenseNDimArray(zeros((dim,)*4).astype(int))
            temp[i, i, k, k] = sp.exp(const * hbar) * coef1[i, i, k, k]
            res *= MatrixTensor2(dim, temp, self.ggs)
        return res

    # I always define the simple roots as \alpha_i = e_{i+1} - e_i
    # The root applied to components not in h is set to be zero. This won't affect anything.
    def root_action_left(self, simple_root_num) -> mat1.MatrixTensor1:
        dim = self.dim
        coef1 = self.coef
        coef = sp.MutableDenseNDimArray(zeros((dim,)*2).astype(int))
        for i in range(dim):
            for k, l in [(x, y) for x in range(dim) for y in range(dim)]:
                if i + 1 == simple_root_num:
                    coef[k, l] -= coef1[i, i, k, l]
                elif i + 1 == simple_root_num % dim + 1:
                    coef[k, l] += coef1[i, i, k, l]
        return mat1.MatrixTensor1(dim, coef)
    
    def root_action_right(self, simple_root_num) -> mat1.MatrixTensor1:
        dim = self.dim
        coef1 = self.coef
        coef = sp.MutableDenseNDimArray(zeros((dim,)*2).astype(int))
        for i in range(dim):
            for k, l in [(x, y) for x in range(dim) for y in range(dim)]:
                if i + 1 == simple_root_num:
                    coef[k, l] -= coef1[k, l, i, i]
                elif i + 1 == simple_root_num % dim + 1:
                    coef[k, l] += coef1[k, l, i, i]
        return mat1.MatrixTensor1(dim, coef)

    # Matrices of GGS type have the form \sum a_{ijk}e_{ij}\otimes e_{k, i+k-j (mod n)}
    def check_of_ggs_type(self):
        dim = self.dim
        coef = self.coef
        for i, j in [(x, y) for x in range(dim) for y in range(dim)]:
            for k, l in [(x, y) for x in range(dim) for y in range(dim)]:
                if l != (i + k - j) % dim and coef[i, j, k, l] != 0:
                    print(f"{i} {j} {k} {l}")
                    return MatrixTensor2(dim, coef, False)
        return MatrixTensor2(dim, coef, True)