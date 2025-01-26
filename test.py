import math, sympy as sp, triple, ybe, mat1
from sympy.combinatorics import Permutation
import mat2, time

# Need to declare the variables first
x = sp.Symbol("x")
h = sp.Symbol("h")

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


para = [5, 4, 0, 7, 8, 0, 0, 0]
trip = triple.BDTriple(para)
s = trip.choose_r0(only_return_s=True)
r0 = s + sp.Rational(1, 2) * mat2.casimir(trip.n)

# ggs_conjecture_aux does NOT use the conjectured value `conjec`.
conjec = sp.Rational(1, 1)
R1 = ybe.ggs_conjecture_aux(trip, x, h, True, conjec)


print(ybe.qybe(R1, x).simplify())
print('=========================================')
print(s)
print('=========================================')
print(R1)


# attention = [0, 0, 1, 2, 0, 3]
# coef1 = ybe.qybe1_aux(R1, x, attention).coef
# coef2 = ybe.qybe2_aux(R1, x, attention).coef
# print(coef1[attention].simplify())
# print(coef[attention].simplify())

coef_print = False

if coef_print:
    coef = ybe.qybe(R1, x).coef
    n = trip.n
    for i, j in [(x, y) for x in range(n) for y in range(n)]:
        for k, l in [(x, y) for x in range(n) for y in range(n)]:
            for p, q in [(x, y) for x in range(n) for y in range(n)]:
                start = time.time()
                temp = coef[i, j, k, l, p, q].simplify()
                if temp != 0:
                    print(f"i={i+1} j={j+1} k={k+1} l={l+1} p={p+1} q={q+1}")
                    print(temp)
                    print(f"The above simplification takes {time.time() - start} seconds")