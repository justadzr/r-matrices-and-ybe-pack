import triple
import sympy as sp
from typing import List, Tuple, Iterable
import nonassoc_affine, ybe, mat2

def cyclic_interval(n: int, start: int, length: int) -> Tuple[int, ...]:
    """
    Consecutive cyclic subset of Z/nZ represented by {1,...,n}.

    start ∈ {1,...,n}
    length ≥ 0

    Returns:
        (start, start+1, ..., start+length-1) mod n,
        with values in {1,...,n}.
    """
    if not (1 <= start <= n):
        raise ValueError("start must be in {1,...,n}")
    if length < 0:
        raise ValueError("length must be nonnegative")

    return tuple(((start - 1 + k) % n) + 1 for k in range(length))


def all_cyclic_intervals(
    n: int,
    include_empty: bool = False,
    include_full: bool = True,
    as_sets: bool = False
) -> List[Iterable[int]]:
    """
    Generate all consecutive subsets (cyclic intervals) of Z/nZ = {1,...,n}.

    Options:
    - include_empty: include the empty interval
    - include_full: include the full set {1,...,n}
    - as_sets: return intervals as sets instead of ordered tuples
    """
    lengths = []
    if include_empty:
        lengths.append(0)
    lengths.extend(range(1, n))
    if include_full:
        lengths.append(n)

    intervals = []
    seen = set()

    for L in lengths:
        for s in range(1, n + 1):
            interval = cyclic_interval(n, s, L)

            if as_sets:
                key = frozenset(interval)
            else:
                # Canonical representative for the full interval
                if L == n:
                    i = interval.index(1)
                    key = interval[i:] + interval[:i]
                else:
                    key = interval

            if key not in seen:
                seen.add(key)
                intervals.append(set(interval) if as_sets else interval)

    return intervals

trip = triple.BDTriple([4, 3, 0, 0])
trip_t = triple.BDTriple([0, 0, 2, 1])

xx = sp.Symbol('x')
qn = sp.Symbol('qn', positive=True)
n = trip.n
qq = qn ** (sp.Rational(n, 2))
twist, twist_s = ybe.ess_twist(trip_t, xx, qn)
formula = ybe.ggs_conjecture_rat(trip, xx, qn)
s = trip.choose_r0(only_return_s=True)
s_t = trip_t.choose_r0(only_return_s=True)

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

to_test = 1 / (qq - qq ** (-1)) * twist_s
temp1 = to_test.subs(xx, 1/xx).subs(qn, 1/qn).transpose()
temp2 = (-s).exp_rat(qn, n, True) * formula * (-s).exp_rat(qn, n, True)
for i, j in [(x, y) for x in range(n) for y in range(n)]:
    for k, l in [(x, y) for x in range(n) for y in range(n)]:
        if (temp1.coef[i, j, k, l] + temp2.coef[i, j, k, l]).simplify() != 0:
            print((i+1, j+1, k+1, l+1))
            print((temp1.coef[i, j, k, l] / temp2.coef[i, j, k, l]).simplify())
            
print(f"Formula twist difference for {trip} ===================================================")
print(( temp1 + temp2  ).simplify())

# print("Rearranged twist=====================================")
# twist_new = (-s).exp_rat(qn, n, True) * (twist) * (-s).exp_rat(qn, n, True)
# twist_new_new = 1/(qq - qq ** (-1)) * twist_new.transpose()

# print(twist_new_new.simplify())
# print("Difference ======================================")

# print((twist_new_new - formula).simplify())