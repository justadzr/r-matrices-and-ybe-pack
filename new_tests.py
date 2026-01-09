import triple
import sympy as sp
from typing import List, Tuple, Iterable
import nonassoc_affine, ybe, mat2, sys


xx = sp.Symbol('x')

qn = sp.Symbol('qn', positive=True)


# _, twist_const = ybe.ess_twist_const(trip, qn)
# formula_const = (-s).exp_rat(qn, n, True) * ybe.ggs_conjecture_constant(trip, qn) * (-s).exp_rat(qn, n, True)
# print((twist_const - formula_const).simplify())
# print("Formula: ===================================================================")
# print(((-s).exp_rat(qn, n, True) * formula * (-s).exp_rat(qn, n, True)).simplify())
# print("Twist: =====================================================================")
# 
# print(to_test.simplify())
# print("===========================================================================")
# print(to_test.subs(xx, 1/xx).subs(qn, 1/qn).transpose().simplify())

for n in range(4, 13):
#     print("*************************************************")
#     print(f"When n = {n}")
#     qq = qn ** (sp.Rational(n, 2))

    # with open(f"codes\\nonassociative-affine-triples-{n}.txt", "r") as f:
    #     unclean_triples = f.read()[2:-2].split('], [')
    #     triples = []
    #     for unclean in unclean_triples:
    #         t = str.maketrans({',': ' ', '[': ' ', ']': ' '})
    #         lst = list(map(int, unclean.translate(t).split()))
    #         if lst:
    #             triples += [triple.BDTriple(lst)]
    for n in range(13, 21):
        triples = nonassoc_affine.nonassoc_affine_triples(n)
        with open(f"nonassociative-affine-triples-{n}.txt", "w") as f1:
            print(triples, file=f1)
            with open(f"server-output/nightmare-for-{n}.txt", "w") as f2:
                for trip in triples:
                    PL, _ = trip.passing_orders_with_orientations()
                    # print(PL)
                    for k in PL:
                        # if k[2] == 1 and PL[k] == 1:
                        #     print(f"For the triple {trip}, the pair {k} is left passed and reversed.")
                        no = True
                        if PL[k] == 1/2:
                            for j in PL:
                                if k[0][1] == j[0][0] and k[1]:
                                    print(f"For the triple {trip}, quadruple {j, k} is a nightmare.", file=f2)
                                    no = False
        # if no:
        #     print(f"{trip} safe.")
    

#     with open(f"codes\\temp_output\\J_factors_{n}.txt", "w") as f:
#         for trip in triples:
#             print("=================================================", file=f)
#             print(f"For the triple {trip}", file=f)
#             J_inv, J21, _, _, std = ybe.ess_twist(trip, xx, qn)
#             print("All J^-1 factors: ------------------------------------------------", file=f)
#             print([(j - mat2.identity(n)).simplify() for j in J_inv], file=f)
#             print("All J21 factors: ------------------------------------------------", file=f)
#             print([(j - mat2.identity(n)).simplify() for j in J21], file=f)

# trip = triple.BDTriple([3, 4, 5, 6, 0, 0])
# focus_lst = [1, 0, 1, 2]
# trip = triple.BDTriple([7, 8, 0, 0, 4, 3, 0, 0, 12, 11, 0, 0])
# focus_lst = [10, 8, 10, 0]
# trip = triple.BDTriple([0, 1, 9, 8, 7, 0, 0, 0, 0])
# focus_lst = list(map(lambda x:x-1, [7, 2, 6, 2]))
# s = trip.choose_r0(only_return_s=True)
# n = trip.n
# print("=================================================")
# print(f"For the triple {trip}")
# _, _, _, res, std = ybe.ess_twist_paper(trip, xx, qn, focus_lst)
# _, _, _, res_2, std = ybe.ess_twist(trip, xx, qn)

# tuple_temp = [0] * n
# for i in range(n):
#     if trip.tuple[i] != 0:
#         tuple_temp[trip.tuple[i] - 1] = i + 1
# trip_t = triple.BDTriple(tuple_temp)
# res_tt = ybe.ggs_conjecture_rat(trip_t, xx, qn)


# print("=================================================")
# print(res.simplify())
# print("=================================================")
# print(res_2.simplify())
# qq = qn ** (sp.Rational(n, 2))
# to_test = 1 / (qq - qq ** (-1)) * res_2
# temp1 = to_test.subs(xx, 1/xx).subs(qn, 1/qn).transpose()
# print("=================================================")
# print(temp1.simplify())
# print("=================================================")
# print(((-s).exp_rat(qn, n, True) * res_tt * (-s).exp_rat(qn, n, True)).simplify())
        
# print(ybe.qybe_rat(s.exp_rat(qn, n, True) * res * s.exp_rat(qn, n, True), xx).simplify())
        
# triples = [triple.BDTriple([0, 0, 2, 1])]
# n = triples[0].n
# qq = qn ** (sp.Rational(n, 2))
# for trip in triples:
#             tuple_temp = [0] * n
#             for i in range(n):
#                 if trip.tuple[i] != 0:
#                     tuple_temp[trip.tuple[i] - 1] = i + 1
#             trip_t = triple.BDTriple(tuple_temp)
#             _, _, _, twist_s, _ = ybe.ess_twist(trip_t, xx, qn)
#             formula = ybe.ggs_conjecture_rat(trip, xx, qn)
#             s = trip.choose_r0(only_return_s=True)
#             s_t = trip_t.choose_r0(only_return_s=True)
#             to_test = 1 / (qq - qq ** (-1)) * twist_s
#             temp1 = to_test.subs(xx, 1/xx).subs(qn, 1/qn).transpose()
#             temp2 = (-s).exp_rat(qn, n, True) * formula * (-s).exp_rat(qn, n, True)
#             print(temp1)
#             res = ( temp1 + temp2  ).simplify()
#             if str(res) != '0':
#                 print(f"Formula twist difference for {trip} FAILED!")
#                 print(( temp1 + temp2  ).simplify())
#             else:
#                 print(f"Formula twist difference for {trip} SUCCESS!")
#                 print(( temp1 + temp2  ).simplify())

# trip = triple.BDTriple([7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 0, 0, 0, 3, 2, 1, 24, 23, 22, 0, 0, 0])

# trip = triple.BDTriple([3, 4, 5, 6, 0, 1, 8, 0])
# print(trip.n)
# print(trip.connected_components())
# PL, PR = trip.passing_orders()
# for i in PR:
#     if PR[i] > 0:
#         print(i, PR[i])