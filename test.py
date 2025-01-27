import math, sympy as sp, triple, ybe, mat1
from sympy.combinatorics import Permutation
import mat2, time, all_triples

# Need to declare the variables first
x = sp.Symbol("x")
h = sp.Symbol("h")

# triples = all_triples.all_triples(4)
# for trip in triples:
para = [0, 0, 0, 3, 2]
trip = triple.BDTriple(para)
s = trip.choose_r0(only_return_s=True)
r0 = s + sp.Rational(1, 2) * mat2.casimir(trip.n)
# ggs_conjecture_aux does NOT use the conjectured value `conjec`.
conjec = sp.Rational(1, 1)
R1 = ybe.ggs_conjecture_aux(trip, x, h, True)
r1 = ybe.to_trigonometric_solution(trip, x, True)
print(R1)
print('=========================================')
print(r1)
print('=========================================')
print(ybe.qybe(R1, x).simplify())

# n = 3
# triple_empty = triple.BDTriple([0] * n)
# std1, std2 = ybe.ggs_conjecture(triple_empty, x, h, True, 0)
# coef_temp = mat2.to_sparray(n, [0] * pow(n, 4))
# for i, j in [(x, y) for x in range(n) for y in range(n)]:
#     if i != j:
#         coef_temp[i, i, j, j] = sp.Rational(1, 2) - sp.Rational((j - i) % n, n)
# s0 = mat2.MatrixTensor2(n, coef_temp, True)
# standard_part = (- s0).exp(h, True) * std1 * (- s0).exp(h, True) + std2
# print(std1)
# print('=========================================')
# print(std2)
# print('=========================================')
# print(standard_part.pr_to_sln())
# print('=========================================')
# print(ybe.to_trigonometric_solution(triple_empty, x, True))

# check a nonzero QYBE(R) by each term
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
                    print(f"The above simplification took {time.time() - start:0.2f} seconds.")