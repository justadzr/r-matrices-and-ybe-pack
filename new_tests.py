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
    print("========================================")
    print(f"When n = {n}")
    qq = qn ** (sp.Rational(n, 2))

    with open(f"nonassociative-affine-triples-{n}.txt", "r") as f:
        unclean_triples = f.read()[2:-2].split('], [')
        triples = []
        for unclean in unclean_triples:
            t = str.maketrans({',': ' ', '[': ' ', ']': ' '})
            lst = list(map(int, unclean.translate(t).split()))
            if lst:
                triples += [triple.BDTriple(lst)]
        for trip in triples:
            trip_t = [0] * n
            for i in range(n):
                if trip[i] != 0:
                    trip_t[trip[i]-1] = i + 1
            twist, twist_s = ybe.ess_twist(trip_t, xx, qn)
            formula = ybe.ggs_conjecture_rat(trip, xx, qn)
            s = trip.choose_r0(only_return_s=True)
            s_t = trip_t.choose_r0(only_return_s=True)
            to_test = 1 / (qq - qq ** (-1)) * twist_s
            temp1 = to_test.subs(xx, 1/xx).subs(qn, 1/qn).transpose()
            temp2 = (-s).exp_rat(qn, n, True) * formula * (-s).exp_rat(qn, n, True)
            res = ( temp1 + temp2  ).simplify()
            if str(res) != '0':
                print(f"Formula twist difference for {trip} FAILED!", file=sys.stderr)
                print(( temp1 + temp2  ).simplify(), file=sys.stderr)
            else:
                print(f"Formula twist difference for {trip} SUCCESS!")
                print(( temp1 + temp2  ).simplify())