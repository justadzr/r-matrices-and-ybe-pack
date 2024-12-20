import mat3, mat2, mat1
import triple as triple
import sympy as sp
from numpy import zeros

# Need to declare the variables first
(u1, u2, u3) = sp.symbols("u1:4")

def subs(num_or_symbol, x, u):
    if isinstance(num_or_symbol, sp.Expr):
        return num_or_symbol.subs(x, u)
    else:
        return num_or_symbol
    
def r_1213(r: mat2.MatrixTensor2, x: sp.Symbol):
    dim = r.dim
    coef = r.coef
    if r.ggs:
        coef_1213 = sp.MutableDenseNDimArray(zeros((dim,)*6).astype(int))
        for i, j in [(x, y) for x in range(dim) for y in range(dim)]:
            for k, l in [(x, y) for x in range(dim) for y in range(dim)]:
                for p in range(dim):
                    q = (i + k + p - j - l) % dim
                    coef_1213[i, j, k, l, p, q] = subs(coef[i, (k + i - l) % dim, k, l], x, u1 - u2) *\
                        subs(coef[(k + i - l) % dim, j, p, q], x, u1 - u3)
        return mat3.MatrixTensor3(dim, coef_1213, True)
    else:
        return r.to_matrixtensor3_12().subs(x, u1 - u2) * r.to_matrixtensor3_13().subs(x, u1 - u3)
    
def r_1312(r: mat2.MatrixTensor2, x: sp.Symbol):
    dim = r.dim
    coef = r.coef
    if r.ggs:
        coef_1312 = sp.MutableDenseNDimArray(zeros((dim,)*6).astype(int))
        for i, j in [(x, y) for x in range(dim) for y in range(dim)]:
            for k, l in [(x, y) for x in range(dim) for y in range(dim)]:
                for p in range(dim):
                    q = (i + k + p - j - l) % dim
                    coef_1312[i, j, k, l, p, q] = subs(coef[i, (l + j - k) % dim, p, q], x, u1 - u3) *\
                        subs(coef[(l + j - k) % dim, j, k, l], x, u1 - u2)
        return mat3.MatrixTensor3(dim, coef_1312, True)
    else:
        return r.to_matrixtensor3_13().subs(x, u1 - u3) * r.to_matrixtensor3_12().subs(x, u1 - u2)
    
def r_1223(r: mat2.MatrixTensor2, x: sp.Symbol):
    dim = r.dim
    coef = r.coef
    if r.ggs:
        coef_1223 = sp.MutableDenseNDimArray(zeros((dim,)*6).astype(int))
        for i, j in [(x, y) for x in range(dim) for y in range(dim)]:
            for k, l in [(x, y) for x in range(dim) for y in range(dim)]:
                for p in range(dim):
                    q = (i + k + p - j - l) % dim
                    coef_1223[i, j, k, l, p, q] = subs(coef[i, j, k, (k + i - j) % dim], x, u1 - u2) *\
                        subs(coef[(k + i - j) % dim, l, p, q], x, u2 - u3)
        return mat3.MatrixTensor3(dim, coef_1223, True)
    else:
        return r.to_matrixtensor3_12().subs(x, u1 - u2) * r.to_matrixtensor3_23().subs(x, u2 - u3)
    
def r_2312(r: mat2.MatrixTensor2, x: sp.Symbol):
    dim = r.dim
    coef = r.coef
    if r.ggs:
        coef_1213 = sp.MutableDenseNDimArray(zeros((dim,)*6).astype(int))
        for i, j in [(x, y) for x in range(dim) for y in range(dim)]:
            for k, l in [(x, y) for x in range(dim) for y in range(dim)]:
                for p in range(dim):
                    q = (i + k + p - j - l) % dim
                    coef_1213[i, j, k, l, p, q] = subs(coef[k, (l + j - i) % dim, p, q], x, u2 - u3) *\
                        subs(coef[i, j, (l + j - i) % dim, l], x, u1 - u2)
        return mat3.MatrixTensor3(dim, coef_1213, True)
    else:
        return r.to_matrixtensor3_23().subs(x, u2 - u3) * r.to_matrixtensor3_12().subs(x, u1 - u2)
    
def r_1323(r: mat2.MatrixTensor2, x: sp.Symbol):
    dim = r.dim
    coef = r.coef
    if r.ggs:
        coef_1323 = sp.MutableDenseNDimArray(zeros((dim,)*6).astype(int))
        for i, j in [(x, y) for x in range(dim) for y in range(dim)]:
            for k, l in [(x, y) for x in range(dim) for y in range(dim)]:
                for p in range(dim):
                    q = (i + k + p - j - l) % dim
                    coef_1323[i, j, k, l, p, q] = subs(coef[i, j, p, (i + p - j) % dim], x, u1 - u3) *\
                        subs(coef[k, l, (i + p - j) % dim, q], x, u2 - u3)
        return mat3.MatrixTensor3(dim, coef_1323, True)
    else:
        return r.to_matrixtensor3_13().subs(x, u1 - u3) * r.to_matrixtensor3_23().subs(x, u2 - u3)
    
def r_2313(r: mat2.MatrixTensor2, x: sp.Symbol):
    dim = r.dim
    coef = r.coef
    if r.ggs:
        coef_2313 = sp.MutableDenseNDimArray(zeros((dim,)*6).astype(int))
        for i, j in [(x, y) for x in range(dim) for y in range(dim)]:
            for k, l in [(x, y) for x in range(dim) for y in range(dim)]:
                for p in range(dim):
                    q = (i + k + p - j - l) % dim
                    coef_2313[i, j, k, l, p, q] = subs(coef[i, j, (k + p - l) % dim, q], x, u1 - u3) *\
                        subs(coef[k, l, p, (k + p - l) % dim], x, u2 - u3)
        return mat3.MatrixTensor3(dim, coef_2313, True)
    else:
        return r.to_matrixtensor3_23().subs(x, u2 - u3) * r.to_matrixtensor3_13().subs(x, u1 - u3)

# Note this function does not simplify the result.
def cybe(r: mat2.MatrixTensor2, x: sp.Symbol) -> mat3.MatrixTensor3:
    return r_1213(r, x) - r_1312(r, x) + r_1223(r, x) - r_2312(r, x) + r_1323(r, x) - r_2313(r, x)


def qybe1(R: mat2.MatrixTensor2, x: sp.Symbol) -> mat3.MatrixTensor3:
    if R.ggs:
        dim = R.dim
        c1 = R.coef
        c = sp.MutableDenseNDimArray(zeros((dim,)*6).astype(int))
        for i, j in [(x, y) for x in range(dim) for y in range(dim)]:
            for k, l in [(x, y) for x in range(dim) for y in range(dim)]:
                for m in range(dim):
                    n = (i + k + m - j - l) % dim
                    for p in range(dim):
                        q = (k + i + m - p - j) % dim
                        c[i, j, k, l, m, n] += subs(c1[i, (k + i - p) % dim, k, p], x, u1 - u2)\
                            * subs(c1[(k + i - p) % dim, j, m, q], x, u1 - u3)\
                            * subs(c1[p, l, q, n], x, u2 - u3)
        return mat3.MatrixTensor3(dim, c, True)
    else:
        R12 = R.subs(x, u1-u2).to_matrixtensor3_12
        R13 = R.subs(x, u1-u3).to_matrixtensor3_13
        R23 = R.subs(x, u2-u3).to_matrixtensor3_23
        return R12 * R13 * R23

def qybe2(R: mat2.MatrixTensor2, x: sp.Symbol) -> mat3.MatrixTensor3:
    if R.ggs:
        dim = R.dim
        c1 = R.coef
        c = sp.MutableDenseNDimArray(zeros((dim,)*6).astype(int))
        for i, j in [(x, y) for x in range(dim) for y in range(dim)]:
            for k, l in [(x, y) for x in range(dim) for y in range(dim)]:
                for m in range(dim):
                    n = (i + k + m - j - l) % dim
                    for p in range(dim):
                        c[i, j, k, l, m, n] += subs(c1[k, p, m, (k + m - p) % dim], x, u2 - u3)\
                            * subs(c1[i, (j + l - p) % dim, (k + m - p) % dim, n], x, u1 - u3)\
                            * subs(c1[(j + l - p) % dim, j, p, l], x, u1 - u2)
        return mat3.MatrixTensor3(dim, c, True)
    else:
        R12 = R.subs(x, u1-u2).to_matrixtensor3_12
        R13 = R.subs(x, u1-u3).to_matrixtensor3_13
        R23 = R.subs(x, u2-u3).to_matrixtensor3_23
        return R23 * R13 * R12

# Note this function does not simplify the result.
def qybe(R: mat3.MatrixTensor3, x: sp.Symbol) -> mat3.MatrixTensor3:
    return qybe1(R, x) - qybe2(R, x)

def check_continuous_datum(trip: triple.BDTriple, r0: mat2.MatrixTensor2):
    dim = r0.dim
    if r0 + r0.swap() - mat2.casimir(trip.n) != mat2.zero(dim):
        return False
    for i in trip.g1:
        if(r0.root_action_left(trip.T(i)) + r0.root_action_right(i) != mat1.zero(dim)):
            print(r0.root_action_left(trip.T(i)) + r0.root_action_right(i))
            return False
    return True

def to_trigonometric_solution(triple: triple.BDTriple, x: sp.Symbol, standard_part: bool) \
    -> mat2.MatrixTensor2:
    n = triple.n
    T = triple.T
    components = triple.connected_components()
    coef1 = sp.MutableDenseNDimArray(zeros((n,)*4).astype(int))
    coef2 = sp.MutableDenseNDimArray(zeros((n,)*4).astype(int))

    # First I choose r0
    r0 = triple.choose_r0(only_return_s = False)

    # Next I insert the standard part (minus its terms in the Cartan)
    for m in range(1, n):
        for i, j in [(x, y) for x in range(n) for y in range(n)]:
            if (i - j) % n == m:
                # The standard part
                coef1[i, j, j, i] += int(standard_part) * \
                    sp.exp(sp.Rational(m, n) * x) / (sp.exp(x) - 1)
                
    def in_one_component(i, j):
        for connected in components:
            p = (i - j) % n
            if set([red(j + k, n) for k in range(p)]).issubset(connected):
                return True
        return False
    
    def take_out_ind(symb):
        s = str(symb)
        if s[0] == "-":
            return int(s[-1]), int(s[2])
        else:
            return int(s[1]), int(s[-1])
    
    def red(a, b):
        return (a - 1) % b + 1

    # At last the nonstandard part
    e = sp.symbols(f"e1:{n + 1}")
    for m in range(1, n):
        for i, j in [(x, y) for x in range(n) for y in range(n)]:
                i_human = i + 1
                j_human = j + 1
                if (i - j) % n == m:
                    while in_one_component(i_human, j_human):
                        a = i_human - 1
                        b = j_human - 1
                        root = 0
                        
                        p = (a - b) % n
                        for q in range(p):
                            root += e[T(red(j_human + q, n)) % n] - e[T(red(j_human + q, n)) - 1]
                        
                        k_human, l_human = take_out_ind(root)
                        k, l = k_human - 1, l_human - 1
                        indicator = triple.C((i + 1, j + 1), (k_human, l_human))
                        if indicator is None:
                            break
                        else:
                            coef2[k, l, j, i] -= (-1) ** (indicator * (abs(a - b) - 1)) * \
                                sp.exp(sp.Rational(m, n) * x)
                            coef2[j, i, k, l] += (-1) ** (indicator * (abs(b - a) - 1)) * \
                                sp.exp(sp.Rational(-m, n) * x)
                            i_human, j_human = k_human, l_human
                    # print(f"i: {i_human}, j: {j_human}, k: {k_human}, l: {l_human}")
                    # print(f"indicator: {indicator}")
                    # print(f"Here counter is {counter}\n======================")

    return  mat2.MatrixTensor2(n, coef2, True) + int(standard_part) * \
        (mat2.MatrixTensor2(n, coef1, True) + (1 / (sp.exp(x) - 1)) * mat2.casimir(n) + r0)

def ggs_conjecture(triple: triple.BDTriple, x: sp.Symbol):
    pass