import math, sympy as sp, triple, ybe, mat1
from sympy.combinatorics import Permutation
import mat2, time, all_triples

# Need to declare the variables first
x = sp.Symbol("x")
h = sp.Symbol("h")
qn = sp.Symbol("qn")

# para = [3, 4, 0, 0]
# trip = triple.BDTriple(para)
# R = ybe.ggs_conjecture_aux_rat(trip, x, q_nth)
# print(ybe.qybe_rat(R, x).simplify())

# Wrong: [2, 0, 0, 5, 6, 1]
# Correct: [0, 5, 0, 2, 3, 0]

trip = triple.BDTriple([2, 0, 0, 5, 6, 1])
n = trip.n
# print(trip.valid())

# # s = trip.choose_r0(only_return_s=True)
# # print(s)

# conj1 = 0
# expo1 = sp.Rational(n, 2) * (conj1 - sp.Rational(19, 14))
# conj2 = 1
# expo2 = sp.Rational(n, 2) * (conj2 - sp.Rational(-11, 7))
# conj3 = 0
# expo3 = sp.Rational(n, 2) * (conj3 - sp.Rational(19, 14))

# R_coef = ybe.ggs_conjecture_rat(trip, x, qn).coef
# R_coef[1, 2, 6, 5] = (1/(qn ** expo1 *x))
# R_coef[6,5, 1, 2] = -(qn ** expo1 *x)
# R_coef[3,4, 6, 5] = (1/(qn ** expo2 *x))
# R_coef[6,5, 3,4] = -(qn ** expo2 *x)
# R_coef[3,4, 1, 0] = (1/(qn ** expo3 *x))
# R_coef[1,0,3,4] = -(qn ** expo3 *x)

# # R_coef[3,5, 6, 4] = (1/(qn ** expo3 *x**2))
# # R_coef[6,4, 3,5] = -(qn ** expo3 *x**2)

# R = mat2.MatrixTensor2(n, R_coef, True)
# print(R)
# QYBER = ybe.qybe_rat(R, x)
# # print(QYBER.simplify())


# for i, j in [(x, y) for x in range(n) for y in range(n)]:
#     for k, l in [(x, y) for x in range(n) for y in range(n)]:
#         for p, q in [(x, y) for x in range(n) for y in range(n)]:
#             temp = QYBER.coef[i,j,k,l,p,q].simplify()
#             if str(temp) != "0":
#                 print(f"i={i} j={j} k={k} l={l} p={p} q={q}:")
#                 print(temp)

# attention = [0, 0, 3, 4, 6, 5]
# ybe.qybe1_rat_aux(R, x, attention)
# ybe.qybe2_rat_aux(R, x, attention)

# triples = all_triples.all_triples(5)
# num = len(triples)
# print([trip.tuple for trip in triples])
# R = ybe.ggs_conjecture_rat(trip, x, qn)
print("===================================================")
R = ybe.ggs_conjecture_rat_new(trip, x, qn)
r = ybe.to_constant_solution(trip, True)
# R_coef = R.coef
# R_coef_new = mat2.to_sparray(n, [0] * pow(n, 4))
# for i, j in [(x, y) for x in range(n) for y in range(n)]:
#     for k, l in [(x, y) for x in range(n) for y in range(n)]:
#         R_coef_new[i, j, k, l] = R_coef[j, i, l, k]
# print("============================")
# R = mat2.MatrixTensor2(n, R_coef_new, True)
# print(R)
# print(ybe.qybe_rat(R, x).simplify_rat())
# # print((R12 * R13 * R23 - R23 * R13 * R12).simplify())
# attention = [0, 0, 0, 2, 0, 2]
# ybe.qybe1_rat_aux(R, x, attention)
# ybe.qybe2_rat_aux(R, x, attention)
# print("============================")
# rb = 1 / (x ** n - 1) * (x **n * r + r.swap())
# print(rb.simplify())
# print("============================")
# RB = R + 1 / (x ** n - 1) * mat2.casimir_gl_not_proj(n)
# print(ybe.qybe_rat(R, x).simplify())
# print("============================")
# conj = ybe.ggs_conjecture_rat(trip, x, qn)
# print(ybe.qybe_rat(conj, x).simplify())
# print("============================")
# RB_coef = RB.coef
# rb_coef = rb.coef
# for i, j in [(x, y) for x in range(n) for y in range(n)]:
#     for k, l in [(x, y) for x in range(n) for y in range(n)]:
#         if i - j == l - k:
#             RB_coef[i, j, k, l] *= x ** ((i - j))
#             rb_coef[i, j, k, l] *= x ** ((i - j))
# # rtsl = ybe.to_trigonometric_solution(trip, x, True) #.pr_to_sln()
# RBsl = mat2.MatrixTensor2(n, RB_coef, True).simplify()
# rbsl = mat2.MatrixTensor2(n, rb_coef, True).simplify()
# print("Generations done.")
# print(rbsl)
# print("============================")
# print(RBsl)

# QYBER = ybe.qybe_rat(RB, x)
# for i, j in [(x, y) for x in range(n) for y in range(n)]:
#     for k, l in [(x, y) for x in range(n) for y in range(n)]:
#         for p, q in [(x, y) for x in range(n) for y in range(n)]:
#             temp = QYBER.coef[i,j,k,l,p,q].simplify()
#             if str(temp) != "0":
#                 print(f"i={i} j={j} k={k} l={l} p={p} q={q}:")
#                 print(temp)


# incorrect = []
# i = 0
# for trip in triples:
#     # if trip.tuple != [6, 5, 0, 2, 1, 0, 0] and trip.tuple != [0, 6, 5, 0, 0, 3, 0]:
#         i += 1
#         print(trip.tuple)
#         R = ybe.ggs_conjecture_rat(trip, x, qn)
#         res = ybe.qybe_rat(R, x).simplify_rat()
#         if str(res) != "0":
#             print(res)
#             incorrect.append(str(trip) + "\n")
#         else:
#             print(f"Checked: {i}/{num}")

# print("The incorrect triples are:")
# print(incorrect)