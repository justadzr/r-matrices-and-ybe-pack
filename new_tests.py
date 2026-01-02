import triple
import sympy as sp
from typing import List, Tuple, Iterable
import nonassoc_affine

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

for n in range(4, 13):
    print("========================================")
    print(f"When n = {n}")
    intervals = all_cyclic_intervals(n, include_empty=False, include_full=False)
    a = sp.symbols([f"α{i}" for i in range(1, n + 1)])

    with open(f"nonassociative-affine-triples-{n}.txt", "r") as f:
        unclean_triples = f.read()[2:-2].split('], [')
        triples = []
        for unclean in unclean_triples:
            t = str.maketrans({',': ' ', '[': ' ', ']': ' '})
            lst = list(map(int, unclean.translate(t).split()))
            if lst:
                triples += [triple.BDTriple(lst)]
    intersection_pairs = []
    a = sp.symbols([f"α{i}" for i in range(1, n + 1)])
    for trip in triples:
        intersection_pairs = []
        printed = False
        g1 = set(trip.g1)
        T = trip.T
        for I in intervals:
            if set(I).issubset(g1):
                expo = None
                temp = I[:]
                for k in range(1, n + 1):
                    if 0 in map(T, temp):
                        break
                    elif n in map(T, temp):
                        target = map(T, temp)
                        expo = k
                        break
                    else:
                        temp = map(T, temp)

                if expo is None:
                    continue
                
                tpp_r = sum([a[i-1] for i in I])
                tpp_l = sum([a[i-1] for i in target])

                temp2 = tpp_l.subs(a[n-1], -sum([a[i] for i in range(0, n-1)]))

                J = []
                for i in range(1, n + 1):
                    if temp2.subs(a[i-1], 0) != temp2:
                        J += [i]

                if set(J).issubset(g1):
                    expo2 = None
                    temp = J[:]
                    for k in range(1, n + 1):
                        if 0 in map(T, temp):
                            break
                        elif n in map(T, temp):
                            target2 = map(T, temp)
                            expo2 = k
                            break
                        else:
                            temp = map(T, temp)

                    if expo2 is None:
                        continue
                    tpm_r = sum([a[i-1] for i in target2])
                    tpm_l = sum([a[i-1] for i in J])

                    intersection_pairs.append((expo, expo2, tpp_l, tpp_r, tpm_l, tpm_r))
        if len(intersection_pairs) > 0:
                print("-----------------------------------------------")
                print(f"For triple: {trip}")
                for xxx in intersection_pairs:
                    expo, expo2, tpp_l, tpp_r, tpm_l, tpm_r = xxx
                    print("TP+ pair:")
                    print(f"k = {expo}")
                    print(f"({tpp_l}, {-tpp_r})")
                    print("TP- pair:")
                    print(f"l = {expo2}")
                    print(f"({-tpm_l}, {tpm_r})")