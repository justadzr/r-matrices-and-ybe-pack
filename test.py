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


para = [3, 4, 0, 0]
# para = [2, 3, 0]
trip = triple.BDTriple(para)
s = trip.choose_r0(only_return_s=True)
R1 = ybe.ggs_conjecture(trip, x, h, True)

print('===============')
# print(mat1.e(4, 1, 1).tensor(mat1.e(4, 3, 3)).pr_to_sln())
print(R1)
print('===============')
# print(ybe.qybe(R1, x).simplify())

# n = trip.n
# coef = ybe.qybe(R1, x).coef
# for i, j in [(x, y) for x in range(n) for y in range(n)]:
#     for k, l in [(x, y) for x in range(n) for y in range(n)]:
#         for p, q in [(x, y) for x in range(n) for y in range(n)]:
#             start = time.time()
#             print(f"i={i+1} j={j+1} k={k+1} l={l+1} p={p+1} q={q+1}")
#             print(coef[i, j, k, l, p, q].simplify())
#             print(f"The above simplification takes {time.time() - start} seconds")