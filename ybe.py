import mat2
import triple
import sympy as sp

# Need to declare the variables first
x = sp.Symbol("x")
(u1, u2, u3) = sp.symbols("u1:4")

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