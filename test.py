import math, sympy as sp, triple, ybe
from sympy.combinatorics import Permutation
from time import time

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


para = [3, 4, 0, 0, 0]
trip = triple.BDTriple(para)
s = trip.choose_r0(only_return_s=True)
print(1-s.root_action(2, 1, 4, 3))
# r = ybe.ggs_conjecture(trip, x, h, small_r=True)
# print(ybe.qybe(r, x).simplify())

# r1 = ybe.to_trigonometric_solution(trip, x, True)
# r2 = ybe.to_trigonometric_solution(trip, x, False)
# print(ybe.cybe(r1, x).simplify())
# print(r2)

# print(r1)
# r2 = mat2.MatrixTensor2(r1.dim, r1.coef, False)
# c = ybe.cybe(r2, x)
# print(c.simplify())
# dim = c.dim
# for i, j in [(x, y) for x in range(dim) for y in range(dim)]:
#     for k, l in [(x, y) for x in range(dim) for y in range(dim)]:
#         for p, q in [(x, y) for x in range(dim) for y in range(dim)]:
#                 print(f"i {i} j {j} k {k} l {l} p {p} q {q}")
#                 print(c.coef[i, j, k, l, p, q].simplify())