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

def check_continuous_datum(trip, r0):
    if str(r0 + r0.swap() - mat2.casimir(trip.n)) != "0":
        return False
    for i in trip.g1:
        if(str(r0.root_action_left(trip.T(i)) + r0.root_action_right(i)) != "0"):
            print(r0.root_action_left(trip.T(i)) + r0.root_action_right(i))
            return False
    return True

# Here r0 is the continuous datum we choose.
# Not including the computation of r0 in this method prevents us from repeatedly calculating r0.
def to_trigonometric_solution(triple, r0, x) -> mat2.MatrixTensor2:
    n = triple.n
    coef = mat2.to_sparray([0] * pow(n, 4))
    for m in range(1, n+1):
        for i, j in [(x, y) for x in range(n) for y in range(n)]:
            if (i - j) % n == m:
                coef[i, j, j, i] += sp.exp(m * x / n) / (sp.exp(x) - 1)
    return mat2.MatrixTensor2(coef) + r0