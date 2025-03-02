import math, sympy as sp, triple, ybe, mat1
from sympy.combinatorics import Permutation
import mat2, time, all_triples

# Need to declare the variables first
x = sp.Symbol("x")
h = sp.Symbol("h")
q_nth = sp.Symbol("qn")

# para = [3, 4, 0, 0]
# trip = triple.BDTriple(para)
# R = ybe.ggs_conjecture_aux_rat(trip, x, q_nth)
# print(ybe.qybe_rat(R, x).simplify())

triples = all_triples.all_triples(7)
num = len(triples)
print(num)
incorrect = []
i = 0
for trip in triples:
    i += 1
    R = ybe.ggs_conjecture_rat(trip, x, q_nth)
    res = ybe.qybe_rat(R, x).simplify_rat()
    if str(res) != "0":
        print(res)
        incorrect.append(str(trip) + "\n")
    else:
        print(f"Checked: {i}/{num}")
print("The incorrect triples are:" + incorrect)

# n =4
# triple_empty = triple.BDTriple([0] * n)
# std1, std2 = ybe.ggs_conjecture(triple_empty, x, h, True, 0)
# coef_temp = mat2.to_sparray(n, [0] * pow(n, 4))
# for i, j in [(x, y) for x in range(n) for y in range(n)]:
#     if i != j:
#         coef_temp[i, i, j, j] = sp.Rational(1, 2) - sp.Rational((j - i) % n, n)
# s0 = mat2.MatrixTensor2(n, coef_temp, True)

# coef_temp1 = mat2.to_sparray(n, [0] * pow(n, 4))
# for k in range(1, n):
#     for i in range(n):
#         coef_temp1[(i + k) % n, (i + k) % n, i, i] = sp.exp(h * (sp.Rational(k, 2 * n) + sp.Rational(1, 2))) / (sp.exp(h) - 1)
# s_test = mat2.MatrixTensor2(n, coef_temp1, True)

# standard_part = (- s0).exp(h, True) * std1 * (- s0).exp(h, True) + std2
# print((- s0).exp(h, True) * std2 * (- s0).exp(h, True))
# print('=========================================')
# print(s_test)
# print('=========================================')
# print(std2)
# print('=========================================')
# print(standard_part)
# print('=========================================')
# print(ybe.to_trigonometric_solution(trip, x, True))

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