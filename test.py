import mat2
import triple
import sympy as sp

# Need to declare the variables first
x = sp.Symbol("x")
(u1, u2, u3) = sp.symbols("u1:4")

# An trigonometric example corresponding to the triple [2 0] of sl(2) in Schedler's notation
coef = [
    (sp.exp(x)+1)/(4*(sp.exp(x)-1)), 0, 0, -(sp.exp(x)+1)/(4*(sp.exp(x)-1)), 
    0, (sp.exp(-x/2)-sp.exp(x/2)), 1/(sp.exp(x/2)-sp.exp(-x/2)), 0, 
    0, 1/(sp.exp(x/2)-sp.exp(-x/2)), 0, 0,
    -(sp.exp(x)+1)/(4*(sp.exp(x)-1)), 0, 0, (sp.exp(x)+1)/(4*(sp.exp(x)-1))
]

# A standard trigonometric solution of sl(3)
# a = 0, b = -1/3
const = sp.Rational(1, 3)
coef1 = [
    const + 2/(3*(sp.exp(x)-1)), 0, 0, 0, - 1/(3*(sp.exp(x)-1)), 0, 0, 0, -const - 1/(3*(sp.exp(x)-1)),
    0, 0, 0, sp.exp(2*const * x)/(sp.exp(x)-1), 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, sp.exp(const * x)/(sp.exp(x)-1), 0, 0,
    0, sp.exp(const * x)/(sp.exp(x)-1), 0, 0, 0, 0, 0, 0, 0,
    -const - 1/(3*(sp.exp(x)-1)), 0, 0, 0, const + 2/(3*(sp.exp(x)-1)), 0, 0, 0, - 1/(3*(sp.exp(x)-1)),
    0, 0, 0, 0, 0, 0, 0, sp.exp(2*const * x)/(sp.exp(x)-1), 0,
    0, 0, sp.exp(2 * const * x)/(sp.exp(x)-1), 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, sp.exp(const * x)/(sp.exp(x)-1), 0, 0, 0,
    - 1/(3*(sp.exp(x)-1)), 0, 0, 0, -const - 1/(3*(sp.exp(x)-1)), 0, 0, 0, const + 2/(3*(sp.exp(x)-1))
]
# Adding the part corresponding to the triple [2 0 0] in Schedler's notation
coef11 = [0] * pow(3, 4)
coef11[mat2.hash_human(3, 3, 2, 1, 2)] = -sp.exp(x/3)
coef11[mat2.hash_human(3, 1, 2, 3, 2)] = sp.exp(-x/3)
r1 = mat2.MatrixTensor2(2, mat2.to_sparray(2, coef))
r2 = mat2.MatrixTensor2(3, mat2.to_sparray(3, coef1))
r3 = r2 + mat2.MatrixTensor2(3, mat2.to_sparray(3, coef11))



# I should write these in the class MatrixTensor2, 
# but in that case we can't substitute explicitly.
def cybe(r):
    r12 = r.subs(x, u1-u2).to_matrixtensor3_12()
    r13 = r.subs(x, u1-u3).to_matrixtensor3_13()
    r23 = r.subs(x, u2-u3).to_matrixtensor3_23()
    return r12.comm(r13) + r12.comm(r23) + r13.comm(r23)

def qybe(R):
    R12 = R.subs(x, u1-u2).to_matrixtensor3_12
    R13 = R.subs(x, u1-u3).to_matrixtensor3_13
    R23 = R.subs(x, u2-u3).to_matrixtensor3_23
    return R12 * R13 * R23 - R23 * R13 * R12

# print(cybe(r3).simplify())

T = triple.BDTriple([5, 4, 0, 0, 0, 0, 6])
print(T)
print(T.connected_components())
print(T.connected_components_img())
print(T.associative())