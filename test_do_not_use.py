import math, sympy as sp, triple, ybe, mat1
from sympy.combinatorics import Permutation
import mat2, time, nonassoc_affine
from anytree import Node

# Need to declare the variables first
x = sp.Symbol("x")
h = sp.Symbol("h")
qn = sp.Symbol("qn")
q = sp.Symbol("q")

# para = [3, 4, 0, 0]
# trip = triple.BDTriple(para)
# R = ybe.ggs_conjecture_aux_rat(trip, x, q_nth)
# print(ybe.qybe_rat(R, x).simplify())

# Wrong: [2, 0, 0, 5, 6, 1]
# Correct: [0, 5, 0, 2, 3, 0]

# trip = triple.BDTriple([7, 8, 9, 10, 11, 12, 0, 0, 0, 0, 0, 0])
# n = trip.n

# # # s = trip.choose_r0(only_return_s=True)
# # # print(s)

# # conj1 = 0
# # expo1 = sp.Rational(n, 2) * (conj1 - sp.Rational(19, 14))
# # conj2 = 1
# # expo2 = sp.Rational(n, 2) * (conj2 - sp.Rational(-11, 7))
# # conj3 = 0
# # expo3 = sp.Rational(n, 2) * (conj3 - sp.Rational(19, 14))

# # R_coef = ybe.ggs_conjecture_rat(trip, x, qn).coef
# # R_coef[1, 2, 6, 5] = (1/(qn ** expo1 *x))
# # R_coef[6,5, 1, 2] = -(qn ** expo1 *x)
# # R_coef[3,4, 6, 5] = (1/(qn ** expo2 *x))
# # R_coef[6,5, 3,4] = -(qn ** expo2 *x)
# # R_coef[3,4, 1, 0] = (1/(qn ** expo3 *x))
# # R_coef[1,0,3,4] = -(qn ** expo3 *x)

# # # R_coef[3,5, 6, 4] = (1/(qn ** expo3 *x**2))
# # # R_coef[6,4, 3,5] = -(qn ** expo3 *x**2)

# # R = mat2.MatrixTensor2(n, R_coef, True)
# # print(R)
# # QYBER = ybe.qybe_rat(R, x)
# # # print(QYBER.simplify())


# # for i, j in [(x, y) for x in range(n) for y in range(n)]:
# #     for k, l in [(x, y) for x in range(n) for y in range(n)]:
# #         for p, q in [(x, y) for x in range(n) for y in range(n)]:
# #             temp = QYBER.coef[i,j,k,l,p,q].simplify()
# #             if str(temp) != "0":
# #                 print(f"i={i} j={j} k={k} l={l} p={p} q={q}:")
# #                 print(temp)

# # attention = [0, 0, 3, 4, 6, 5]
# # ybe.qybe1_rat_aux(R, x, attention)
# # ybe.qybe2_rat_aux(R, x, attention)
# # nonassoc_affine.nonassoc_affine_triples(5)
# # triples = all_triples.all_triples(5)
# # num = len(triples)
# # # print([trip.tuple for trip in triples])
# # R = ybe.ggs_conjecture_rat_new(trip, x, qn)
# # # print("===================================================")
# R = ybe.ggs_conjecture_constant(trip, qn)
# # r = ybe.to_constant_solution(trip, True)
# # # R_coef = R.coef
# # # R_coef_new = mat2.to_sparray(n, [0] * pow(n, 4))
# # # for i, j in [(x, y) for x in range(n) for y in range(n)]:
# # #     for k, l in [(x, y) for x in range(n) for y in range(n)]:
# # #         R_coef_new[i, j, k, l] = R_coef[j, i, l, k]
# # # print("============================")
# # # R = mat2.MatrixTensor2(n, R_coef_new, True)
# # # print(R)
# # print(ybe.qybe_rat(R, x).simplify_rat())
# # # print((R12 * R13 * R23 - R23 * R13 * R12).simplify())
# # attention = [0, 0, 0, 2, 0, 2]
# # ybe.qybe1_rat_aux(R, x, attention)
# # ybe.qybe2_rat_aux(R, x, attention)
# # print("============================")
# # rb = 1 / (x ** n - 1) * (x **n * r + r.swap())
# # # print(rb.simplify())
# # print("============================")
# # RB = R + 1 / (x ** n - 1) * mat2.casimir_gl_not_proj(n)
# # print(ybe.qybe_rat(R, x).simplify())
# # # print("============================")
# # # conj = ybe.ggs_conjecture_rat(trip, x, qn)
# # # print(ybe.qybe_rat(conj, x).simplify())
# # # print("============================")
# # RB_coef = RB.coef
# # # # rb_coef = rb.coef
# # for i, j in [(x, y) for x in range(n) for y in range(n)]:
# #     for k, l in [(x, y) for x in range(n) for y in range(n)]:
# #         if i - j == l - k:
# #             RB_coef[i, j, k, l] *= x ** ((i - j))
# # #             # rb_coef[i, j, k, l] *= x ** ((i - j))
# # # # rtsl = ybe.to_trigonometric_solution(trip, x, True) #.pr_to_sln()
# # RB = mat2.MatrixTensor2(n, RB_coef, True).simplify()
# # R = ybe.ggs_conjecture_rat_new(trip, x, qn)
# # # # rbsl = mat2.MatrixTensor2(n, rb_coef, True).simplify()
# # # # print("Generations done.")
# # # # print(rbsl)
# # print(RB)
# # print("============================")
# # print(R)
# # print("============================")
# # print((RB - R).simplify())


# # incorrect = []
# # i = 0
# # for trip in triples:
# #     # if trip.tuple != [6, 5, 0, 2, 1, 0, 0] and trip.tuple != [0, 6, 5, 0, 0, 3, 0]:
# #         i += 1
# #         print(trip.tuple)
# #         R = ybe.ggs_conjecture_rat(trip, x, qn)
# #         res = ybe.qybe_rat(R, x).simplify_rat()
# #         if str(res) != "0":
# #             print(res)
# #             incorrect.append(str(trip) + "\n")
# #         else:
# #             print(f"Checked: {i}/{num}")

# # print("The incorrect triples are:")
# # print(incorrect)

# # trip = triple.BDTriple([4, 3, 0, 6, 1, 0])
# # n = trip.n
# R = ybe.ggs_conjecture_rat(trip, x, qn)
# print("Here")
# QYBER = ybe.qybe_rat(R, x)
# print("simplify starts")
# for i, j in [(x, y) for x in range(n) for y in range(n)]:
#     for k, l in [(x, y) for x in range(n) for y in range(n)]:
#         for p, q in [(x, y) for x in range(n) for y in range(n)]:
#             t = time.time()
#             # print(f"Before: {t}")
#             temp = QYBER.coef[i,j,k,l,p,q].ratsimp()
#             # print(f"After: {time.time()-t}")
#             if str(temp) != "0":
#                 print(f"i={i} j={j} k={k} l={l} p={p} q={q}:")
#                 print(temp)

# R = mat2.zero(5)
# coef = R.coef

# for i in range(5):
#     coef[i,i,i,i]=q

# for (i, j) in [(x,y) for x in range(5) for y in range(5)]:
#     if i != j:
#         coef[i,i,j,j]=1
#     if i > j:
#         coef[i,j,j,i]=q-q ** (-1)
# Rst = mat2.MatrixTensor2(5, coef, True)

# J21_coef = mat2.identity(5).coef

# J21_coef[1,0,3,4]+=q-q ** (-1)
# J21_coef[2,1,2,3]+=q ** sp.Rational(1,2) * (q-q ** (-1))
# J21_coef[2,0,2,4]+=q ** sp.Rational(1,2) * (-q) ** (-1) * (q-q ** (-1))
# J21 = mat2.MatrixTensor2(5, J21_coef, False)
# N = J21.swap()-mat2.identity(5)
# Jinv = mat2.identity(5) - N + N * N - N * N * N

# Ji_coef = mat2.identity(5).coef
# Ji_coef[3,4,1,0] -= (q-q ** (-1))
# Ji_coef[2,3,2,1] -= q ** sp.Rational(1,2) * (q-q ** (-1))
# Ji_coef[2,4,2,0] -= q ** sp.Rational(1,2) * (-q) * (q-q ** (-1))
# Ji = mat2.MatrixTensor2(5,Ji_coef,False)
# print((Ji * J21.swap()).simplify())
# print()

# print((Jinv * Rst * J21).simplify())
# print()

# trip = triple.BDTriple([4, 3, 0, 0, 0])
# n=trip.n
# trip_const = triple.BDTriple([4, 3, 0, 0, 0])
# R = ybe.ggs_conjecture_rat(trip, x, qn)
# R_const = ybe.ggs_conjecture_constant(trip_const, qn)
# s = trip.choose_r0(True)
# R_twist = (-s).exp_rat(qn, n, True) * R_const * (-s).exp_rat(qn, n, True)
# print((((qn**5 - 1) / qn ** sp.Rational(5, 2)) * R_twist).simplify())
# RB = R_const + 1 / (x ** n - 1) * mat2.casimir_gl_not_proj(n)
# RB_coef = RB.coef
# # # rb_coef = rb.coef
# for i, j in [(x, y) for x in range(n) for y in range(n)]:
#     for k, l in [(x, y) for x in range(n) for y in range(n)]:
#         if i - j == l - k:
#             RB_coef[i, j, k, l] *= x ** ((i - j))
# #             # rb_coef[i, j, k, l] *= x ** ((i - j))
# # # rtsl = ybe.to_trigonometric_solution(trip, x, True) #.pr_to_sln()
# RB = mat2.MatrixTensor2(n, RB_coef, True).simplify()
# print(R)

trip = triple.BDTriple([0, 0, 0, 0, 0, 0, 0, 0])
trip_const = triple.BDTriple([0, 0, 0, 0, 0, 0, 0, 0])
n = trip.n
# print(ybe.ggs_conjecture_rat(trip, x, qn))
print("-------------------------")
# print(ybe.ggs_conjecture_constant(trip, qn))
R = ybe.ggs_conjecture_rat_t(trip, x, qn)
s1 = trip.choose_r0(only_return_s=True)
s2 = trip_const.choose_r0(only_return_s=True)
R_const = ybe.ggs_conjecture_constant_sch_s(trip_const, qn ** (-1), s2)
print(s1)
print("-------------------------")
print(s2)

RB = (-1) * (R_const + x ** n / (1 - x ** n) * mat2.casimir_gl_not_proj(n))

RB_coef = RB.coef
for i, j in [(x, y) for x in range(n) for y in range(n)]:
    for k, l in [(x, y) for x in range(n) for y in range(n)]:
        RB_coef[i, j, k, l] *= x ** (j - i)
RB = mat2.MatrixTensor2(n, RB_coef, True).simplify()

print("-------------------------")
print(R)
print("-------------------------")
print(RB.simplify())
print("-------------------------")
print(ybe.qybe_rat(R, x).simplify())
print("-------------------------")
print(ybe.qybe_rat(RB, x).simplify())
print("-------------------------")
print((R - RB).simplify())


R1 = 1 / (1 / (qn ** sp.Rational(n, 2) - qn ** (sp.Rational(-n, 2))) + 1 / (sp.exp(sp.Rational(n, 2) * x) - sp.exp(sp.Rational(-n, 2) * x))) * R
R2 = R1.subs(x, -x).swap()
# print(R1)
# print("-------------------------")
# print(R2)
# print("-------------------------")
# print((R1*R2).simplify())
# print("-------------------------")
# print((R * RB).simplify())
# print("-------------------------")
# print(ybe.qybe_rat(R_const, x).simplify())


def red(a, b):
        return (a - 1) % b + 1