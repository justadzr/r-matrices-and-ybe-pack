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

trip = triple.BDTriple([4, 3, 0, 0])
trip_f = triple.BDTriple([4, 3, 0, 0, 0])
n = trip.n
# print(trip.connected_components())
# print(Y_runs(trip))

R = ybe.ggs_conjecture_rat(trip, sp.exp(sp.Rational(1, n) * x), qn)
print(R)
print("=========================================")

R_norm = (1/((qn ** sp.Rational(n, 2) - qn ** sp.Rational(-n, 2)) ** (-1) + (sp.exp(sp.Rational(1, 2) * x) - sp.exp(sp.Rational(-1, 2) * x)) ** (-1))) * R
print((  R_norm * (  R_norm.subs(x, -x).swap()  )  ).simplify())
print("=========================================")

coeff_ind_in_x_series = -1

R_const = R.expand(x, 0, 10).extract_coeff(x, coeff_ind_in_x_series).subs(q, sp.exp(h / (2 * n)))
print(f"The coefficient of x^{coeff_ind_in_x_series}:")
print(R_const)

# print("=========================================")
# print("r:")
# r = R_const.expand(h, 0, 10).extract_coeff(h, 1)
# print(r)
# print("=========================================")
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