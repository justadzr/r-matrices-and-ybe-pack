import math, sympy as sp, triple, ybe, mat1
from sympy.combinatorics import Permutation
import mat2, time, nonassoc_affine
from anytree import Node
from gsv_basis import *

# Need to declare the variables first
x = sp.Symbol("x")
h = sp.Symbol("h")
qn = sp.Symbol("qn")
q = sp.Symbol("q")
λ = sp.symbols('λ')

trip = triple.BDTriple([0, 0, 0, 0])
# trip_f = triple.BDTriple([4, 3, 0, 0, 0])
n = trip.n
# print(trip.connected_components())
# print(Y_runs(trip))

R = ybe.ggs_conjecture_rat(trip, sp.exp(x/n), qn)
R_match = R.subs(x, -x).transpose().simplify()
qqq = mat2.casimir_gl(n).exp(h/2, True)
qqq_inv = mat2.casimir_gl(n).exp(-h/2, True)

print(R_match)
print("=========================================")
print(qqq)
print("=========================================")
print(qqq_inv)
print("=========================================")
print((qqq_inv * R_match * qqq).simplify() - R_match)
print("=========================================")
print((qqq_inv * (mat1.e(4, 1, 2).tensor(mat1.e(4, 2, 1))) * qqq).simplify())



# R_norm = (1/((qn ** sp.Rational(n, 2) - qn ** sp.Rational(-n, 2)) ** (-1) + (sp.exp(sp.Rational(1, 2) * x) - sp.exp(sp.Rational(-1, 2) * x)) ** (-1))) * R
# print((  R_norm * (  R_norm.subs(x, -x).swap()  )  ).simplify())
print("=========================================")

# coeff_ind_in_x_series = 0

# # R_const = R.expand(x, 0, 10).extract_coeff(x, coeff_ind_in_x_series).subs(qn, sp.exp(sp.Rational(1, n) * h))

# R_const = R.subs(qn, sp.exp(sp.Rational(1, n) * h))

# # print(f"Expanding R using wrt x")
# # print(R_const.expand(x, 0, 4))
# print("Limit at x to infinity")
# print(R_const.limit(x, sp.oo))

# expansion_variable = x
# expansion_terms = 5

# R_const_expansion = R_const.expand(expansion_variable, 0, expansion_terms + 1)

# for i in range(-expansion_terms, expansion_terms + 1):
#         print("=========================================")
#         print(f"at {expansion_variable}^{i}")
#         Ri = R_const_expansion.extract_coeff(expansion_variable, i)
#         print(Ri)
#         print("check cybe:")
#         print(ybe.cybe(Ri, x).simplify())
#         print("=========================================")
# print("r + r^21:")
# print(r + r.swap())

# print("=========================================")

# print("CYBE of r:")
# print(ybe.cybe(r, x))

# print("=========================================")
# print(f"QYBE of the coefficient of x^{coeff_ind_in_x_series}:")

# print(ybe.qybe(R_const, x).simplify())

# print("=========================================")
# rx = ybe.to_trigonometric_solution(trip, x, True)
# print((rx + rx.subs(x, 1/x)).simplify())