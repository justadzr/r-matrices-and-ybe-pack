import triple
import sympy as sp
from typing import List, Tuple, Iterable
import nonassoc_affine, mat1, ybe

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

def red(n, a):
    return ((a - 1) % n) + 1

def e_root(n, alpha, sgn):
    temp = sorted([x-1 for x in alpha])
    res = mat1.e(n, temp[0], (temp[0] + 1) % n)    
    for i in range(1, len(temp)):
        res *= mat1.e(n, temp[i], (temp[i] + 1) % n)
    
    if sgn == 1:
        return res
    else:
        return res.transpose()

xx = sp.Symbol('x')
qn = sp.Symbol('qn')

def to_root_human(n, alpha, a):
    res = 0
    for i in alpha:
        res += a[i-1]
    return res



def Phi(n, conv):
    sgn_x = conv["x sign"]
    x_s = conv["x source"]
    x_t = conv["x target"]
    y_s = conv["y source"]
    y_t = conv["y target"]

    len_x = abs(x_s[0] - x_s[1])
    len_y = abs(y_s[0] - y_s[1])
    
    if len_x < len_y:
        if len_x % (len_y - len_x) == 0:
            j = x_s[1]
            k = x_t[1]
            u_s = (k, red(n, k+len_x-len_y))
            u_t = (red(n, j+len_y), red(n, j+len_x))

            # return({
            #     "x sign": -sgn_x,
            #     "x source": u_s,
            #     "x target": u_t,
            #     "y source": v_s,
            #     "y target": v_t
            # })
            print("Can't deal with BCTQo^i at the moment!")
            return 0
        else:
            j = x_s[1]
            k = x_t[1]
            q = len_x % (len_y - len_x)
            u_s = (k, red(n, k+len_x-len_y+q))
            u_t = (red(n, j+len_y-q), red(n, j+len_x))
            v_s = (red(n, k+len_x-len_y+q), red(n, k+len_x-len_y))
            v_t = (red(n, j+len_y), red(n, j+len_y-q))

            return({
                "x sign": -sgn_x,
                "x source": u_s,
                "x target": u_t,
                "y source": v_s,
                "y target": v_t
            })

    if len_x > len_y:
        if len_x % (len_x - len_y) == 0:
            print("Can't deal with BCTQo^d at the moment!")
            return 0
        else:
            j = x_s[1]
            k = x_t[1]
            q = len_x % (len_y - len_x)
            u_s = (red(n, j+len_x-q), red(n, j+len_y))
            u_t = (red(n, j+len_x-len_y), red(n, k+q))
            v_t = (red(n, k+q), k)
            v_s = (red(n, j+len_y), red(n, j+len_x-q))

            return({
                "x sign": sgn_x,
                "x source": v_s,
                "x target": v_t,
                "y source": u_s,
                "y target": u_t
            })


print("Job starts...")

for n in range(8, 13):
    # n = 11
    print("========================================")
    print(f"When n = {n}")
    intervals = all_cyclic_intervals(n, include_empty=False, include_full=False)

    with open(f"codes\\nonassociative-affine-triples-{n}.txt", "r") as f:
        unclean_triples = f.read()[2:-2].split('], [')
        triples = []
        for unclean in unclean_triples:
            t = str.maketrans({',': ' ', '[': ' ', ']': ' '})
            lst = list(map(int, unclean.translate(t).split()))
            if lst:
                triples += [triple.BDTriple(lst)]

    # triples = [triple.BDTriple([5, 6, 7, 0, 9, 10, 11, 0, 0, 0, 0])]
    intersection_pairs = []
    a = sp.symbols([f"α{i}" for i in range(1, n + 1)])
    for trip in triples:
        T_pairs = []
        neutral = []
        printed = False
        g1 = set(trip.g1)
        T = trip.T
        for I in intervals:
            if set(I).issubset(g1):
                neutral += [(list(I), list(I), 0)]
                expo = None
                temp = I[:]
                for k in range(1, n + 1):
                    if 0 in list(map(T, temp)):
                        break
                    else:
                        target = list(map(T, temp))
                        T_pairs.append((list(I), target, k))
                        temp = list(map(T, temp))

                tpp_r = sum([a[i-1] for i in I])
                tpp_l = sum([a[i-1] for i in target])
        
        PTP = [(y, x, z) for (x, y, z) in T_pairs]
        NTP = neutral + T_pairs[:]

        CTQo = []
        CTQo_a_aux = []
        CTQs = []
        CTQo_a = []
        CTQo_machine = []
        Conv = []

        dic = {1: PTP, -1: NTP}
        signs = [(1, 1), (1, -1), (-1, 1), (-1, -1)]

        for x, y in signs:
            p1 = dic[x]
            p2 = dic[y]
            for alpha, beta, k in p1:
                for gamma, delta, k in p2:
                    e1 = e_root(n, alpha, x)
                    e2 = e_root(n, beta, -x)
                    e3 = e_root(n, gamma, y)
                    e4 = e_root(n, delta, -y)

                    if str(e1 * e3) != '0' and str(e2 * e4) != '0':
                        if x == y:
                            CTQs.append(((x * to_root_human(n, alpha, a), (-x) * to_root_human(n, beta, a)), (y * to_root_human(n, gamma, a), (-y) * to_root_human(n, delta, a))))
                        else:
                            CTQo.append(((x * to_root_human(n, alpha, a), (-x) * to_root_human(n, beta, a)), (y * to_root_human(n, gamma, a), (-y) * to_root_human(n, delta, a))))
                            if x > 0:
                                Conv.append({"x sign": x, "x source": (red(n, beta[-1]+2-1), beta[0]), "x target": (red(n, alpha[-1]+2-1), alpha[0]), "y source": (red(n, gamma[-1]+2-1), gamma[0]), "y target": (red(n, delta[-1]+2-1), delta[0])})

                                CTQo_machine.append(((red(n, beta[-1]+2-1), beta[0]), (red(n, alpha[-1]+2-1), alpha[0]), (red(n, gamma[-1]+2-1), gamma[0]), (red(n, delta[-1]+2-1), delta[0])))
                            else:
                                Conv.append({"x sign": x, "x source": (red(n, alpha[-1]+2-1), alpha[0]), "x target": (red(n, beta[-1]+2-1), beta[0]), "y source": (red(n, delta[-1]+2-1), delta[0]), "y target": (red(n, gamma[-1]+2-1), gamma[0])})

                                CTQo_machine.append(((red(n, alpha[-1]+2-1), alpha[0]), (red(n, beta[-1]+2-1), beta[0]), (red(n, delta[-1]+2-1), delta[0]), (red(n, gamma[-1]+2-1), gamma[0])))

                            if len(alpha) == len(gamma):
                                CTQo_a.append(((x * to_root_human(n, alpha, a), (-x) * to_root_human(n, beta, a)), (y * to_root_human(n, gamma, a), (-y) * to_root_human(n, delta, a))))
                                CTQo_a_aux.append((e1.tensor(e2), e3.tensor(e4)))

        if len(CTQo) > 0:
            print(f"For the triple {trip} we have: ")
            print("Same signs compatible pairs: ")
            print(CTQs)
            print("Opposite signs compatible pairs: ")
            print(CTQo)
            print("Opposite pairs with the same lengths: ")
            print(CTQo_a)
            print("Opposite matrices with the same lengths: ")
            print(CTQo_a_aux)
            print("--------------------------------------------------------")
            PL, PR = trip.passing_orders()
            # dic, _ = ybe.ggs_conjecture_rat_passing_ord(trip, xx, qn)

            def po_get(PO, key):
                if key[0] != key[1]:
                    return PO[(key[0], key[1])]
                else:
                    return 0

            for conv in Conv:
                phi = Phi(n, conv)
                if phi != 0:
                    print(f"(x, y): {conv}")
                    print(f"(u, v): {phi}")
                    plv = po_get(PL, (phi["y source"], phi["y target"]))
                    prv = po_get(PR, (phi["y source"], phi["y target"]))
                    plu = po_get(PL, (phi["x source"], phi["x target"]))
                    pru = po_get(PR, (phi["x source"], phi["x target"]))
                    plx = po_get(PL, (conv["x source"], conv["x target"]))
                    prx = po_get(PR, (conv["x source"], conv["x target"]))
                    ply = po_get(PL, (conv["y source"], conv["y target"]))
                    pry = po_get(PR, (conv["y source"], conv["y target"]))
                    print(f"Pr(v): {prv} Pl(v): {plv}")
                    print(f"Pr(u): {pru} Pl(u): {plu}")
                    print(f"Pr(y): {pry}; Pl(y): {ply}")
                    print(f"Pr(x): {prx}; Pl(x): {plx}")
                    
                    if plv == pry + plx and ply == 0:
                        print("Case (a)")
                    elif plv == ply and plx == 0 and pry == 0:
                        print("Case (b)")
                    else:
                        print("FAILED!")

                    if prv == prx + pru and plu == 0:
                        print("Case (c)")
                    elif prv == plu and pru == 0 and prx == 0:
                        print("Case (d)")
                    else:
                        print("FAILED!")
                    


            # for tp11, tp12, tp21, tp22 in CTQo_machine:
            #     if tp11 != tp12 and dic[(tp11, tp12)] > 1:
            #         print(f"The passing order at {(tp11, tp12)} in the compatible quadruple {((tp11, tp12), (tp21, tp22))} is > 1")
            #     if tp21 != tp22 and dic[(tp21, tp22)] > 1:
            #         print(f"The passing order at {(tp21, tp22)} in the compatible quadruple {((tp11, tp12), (tp21, tp22))} is > 1")
            