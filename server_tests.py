from __future__ import annotations
import os
import sys, mat2
import sympy as sp
from concurrent.futures import ProcessPoolExecutor, as_completed


from typing import List, Tuple

Root = Tuple[int, int]

import triple
import ybe

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

xx = sp.Symbol("x", positive=True)
qn = sp.Symbol("qn", positive=True)
h = sp.Symbol("h", positive=True)


def parse_triples_file(path: str):
    with open(path, "r") as f:
        raw = f.read()
    unclean_triples = raw[2:-2].split('], [')
    out = []
    for unclean in unclean_triples:
        t = str.maketrans({',': ' ', '[': ' ', ']': ' '})
        lst = list(map(int, unclean.translate(t).split()))
        if lst:
            out.append(lst)
    return out


def worker(n: int, lst):
    qq = qn ** (sp.Rational(n, 2))

    trip = triple.BDTriple(lst)

    tuple_temp = [0] * n
    for i in range(n):
        if trip.tuple[i] != 0:
            tuple_temp[trip.tuple[i] - 1] = i + 1
    trip_t = triple.BDTriple(tuple_temp)

    _, _, _, twist_s = ybe.ess_twist(trip_t, xx, qn)
    formula = ybe.ggs_conjecture_rat(trip, xx, qn)

    s = trip.choose_r0(only_return_s=True)

    to_test = twist_s
    temp1 = 1 / (qq  ** (-1) - qq)* (to_test.subs(xx, 1/xx).subs(qn, 1/qn).transpose())
    temp2 = (-s).exp_rat(qn, n, True) * (formula) * (-s).exp_rat(qn, n, True)

    res = (temp1 + temp2).simplify()
    ok = (str(res) == "0")
    # print("==============================================")
    # print("YBE check for the unmodified twist")
    # # print(ybe.qybe_rat((s).exp_rat(qn, n, True) * twist * (s).exp_rat(qn, n, True), xx).simplify())
    # # # print("YBE check for the modified twist")
    # # # print(ybe.qybe_rat((s).exp_rat(qn, n, True) * temp1 * (s).exp_rat(qn, n, True), xx).simplify())
    # print(f"For the inverse triple: {trip_t}, the twist is\n{temp1.simplify()}")
    # print("==============================================")
    # print(f"For the triple: {trip}, the formula is\n{temp2.simplify()}")
    # print("==============================================")
    # u = sp.Symbol("u")
    # ex = xx ** sp.Rational(n)
    # # formula_h_part = (1 / (1/(qq-qq**(-1)) + 1/(xx ** sp.Rational(n, 2)-xx ** sp.Rational(-n, 2))) *formula).pr_to_sln().subs(qn, sp.exp(h/n)).expand(h, 0, 4).extract_coeff(h, 1)
    # # print(formula_h_part.simplify())
    # print("==============================================")
    
    # print(ybe.to_trigonometric_solution(trip, xx, True))
    # print("Focus:")
    # car = mat2.identity(n)
    # for j in J_inv[::-1]:
    #     car = car * j
    #     print((1 / (qq  ** (-1) - qq)*car).simplify())
    if ok:
        print(f"SUCCESS for {trip}")
    else:
        print(f"FAILED for {trip}")
    return str(trip), ok, str(res) 

def worker_const(n: int, lst):
    qq = qn ** (sp.Rational(n, 2))

    trip = triple.BDTriple(lst)

    tuple_temp = [0] * n
    for i in range(n):
        if trip.tuple[i] != 0:
            tuple_temp[trip.tuple[i] - 1] = i + 1
    trip_t = triple.BDTriple(tuple_temp)

    _, twist_s = ybe.ess_twist_const(trip, qn)
    
    formula = ybe.ggs_conjecture_constant(trip, qn)

    s = trip.choose_r0(only_return_s=True)

    to_test = twist_s
    temp1 = to_test
    temp2 = (-s).exp_rat(qn, n, True) * formula * (-s).exp_rat(qn, n, True)
    print("==============================================")
    print(f"For the triple: {trip}, the formula is\n{temp2}")
    print("==============================================")
    print(f"For the inverse triple: {trip}, the modified twist is\n{temp1}")
    print("==============================================")

    res = (temp1 - temp2).simplify()
    ok = (str(res) == "0")
    print("Difference")
    print(res)
    print("==============================================")

# print(triple.BDTriple([4, 3, 1, 0]).passing_orders())

trip_lst = [7, 8, 9, 10, 11, 12, 13, 14, 15, 16,17,18,0,0,0,3, 2, 1, 0]
# print(triple.BDTriple(trip_lst).valid())
worker(len(trip_lst), trip_lst)


# for n in range(4, 9):
#     triples_as_lists = parse_triples_file(f"codes\\nonassociative-affine-triples-{n}.txt")
#     for trip in triples_as_lists:
#         print("=========================================")
#         worker(n, trip)


# def run_all(n_min=4, n_max=12, max_workers=50):
#     workers = min(max_workers, os.cpu_count() or max_workers)

#     for n in range(n_min, n_max + 1):
#         print("========================================")
#         print(f"When n = {n}")

#         triples_as_lists = parse_triples_file(f"nonassociative-affine-triples-{n}.txt")
#         success_path = f"server-output/success-{n}.out"
#         failed_path = f"server-output/failed-{n}.err"

#         with open(success_path, "w") as f_ok, open(failed_path, "w") as f_bad:
#             with ProcessPoolExecutor(max_workers=workers) as ex:
#                 futures = [ex.submit(worker, n, lst) for lst in triples_as_lists]

#                 for fut in as_completed(futures):
#                     trip_str, ok, res_str = fut.result()
#                     if ok:
#                         msg = f"Formula twist difference for {trip_str} SUCCESS!\n{res_str}\n"
#                         print(msg, end="")
#                         f_ok.write(msg + "\n")
#                     else:
#                         msg = f"Formula twist difference for {trip_str} FAILED!\n{res_str}\n"
#                         print(msg, end="", file=sys.stderr)
#                         f_bad.write(msg + "\n")


# if __name__ == "__main__":
#     run_all()