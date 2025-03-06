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

triples = all_triples.nonassociative_triples(11)
num = len(triples)
print(num)

# incorrect = []
# i = 0
# for trip in triples:
#     i += 1
#     R = ybe.ggs_conjecture_rat(trip, x, q_nth)
#     res = ybe.qybe_rat(R, x).simplify_rat()
#     if str(res) != "0":
#         print(res)
#         incorrect.append(str(trip) + "\n")
#     else:
#         print(f"Checked: {i}/{num}")

# print("The incorrect triples are:")
# print(incorrect)