import mat2
import triple
import sympy as sp
import ybe

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

# print(cybe(r3).simplify())
para = [2, 0, 1]
trip = triple.BDTriple(para)
r0 = trip.choose_r0()
casimir = mat2.casimir(len(para))
print(r0)
print(ybe.check_continuous_datum(trip, r0))