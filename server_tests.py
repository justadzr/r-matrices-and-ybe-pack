import os
import sys
import sympy as sp
from concurrent.futures import ProcessPoolExecutor, as_completed

import triple
import ybe

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

xx = sp.Symbol("x")
qn = sp.Symbol("q")   


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

    trip_t = [0] * n
    for i in range(n):
        if trip[i] != 0:
            trip_t[trip[i] - 1] = i + 1

    _, twist_s = ybe.ess_twist(trip_t, xx, qn)
    formula = ybe.ggs_conjecture_rat(trip, xx, qn)

    s = trip.choose_r0(only_return_s=True)

    to_test = twist_s / (qq - qq**(-1))
    temp1 = to_test.subs(xx, 1/xx).subs(qn, 1/qn).transpose()
    temp2 = (-s).exp_rat(qn, n, True) * formula * (-s).exp_rat(qn, n, True)

    res = (temp1 + temp2).simplify()
    ok = (str(res) == "0")

    return n, str(trip), ok, str(res)


def run_all(n_min=4, n_max=12, max_workers=90):
    workers = min(max_workers, os.cpu_count() or max_workers)

    for n in range(n_min, n_max + 1):
        print("========================================")
        print(f"When n = {n}")

        triples_as_lists = parse_triples_file(f"nonassociative-affine-triples-{n}.txt")
        success_path = f"success-{n}.out"
        failed_path = f"failed-{n}.err"

        with open(success_path, "w") as f_ok, open(failed_path, "w") as f_bad:
            with ProcessPoolExecutor(max_workers=workers) as ex:
                futures = [ex.submit(worker, n, lst) for lst in triples_as_lists]

                for fut in as_completed(futures):
                    _, trip_str, ok, res_str = fut.result()
                    if ok:
                        msg = f"Formula twist difference for {trip_str} SUCCESS!\n{res_str}\n"
                        print(msg, end="")
                        f_ok.write(msg + "\n")
                    else:
                        msg = f"Formula twist difference for {trip_str} FAILED!\n{res_str}\n"
                        print(msg, end="", file=sys.stderr)
                        f_bad.write(msg + "\n")


if __name__ == "__main__":
    run_all()