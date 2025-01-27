import random, triple
from itertools import combinations

def all_triples(n : int) -> list[triple.BDTriple]:
    res = []
    lst = [(x + 1) for x in range(n)]
    for i in range(n):
        sublists = []
        sublists.extend(combinations(lst, i + 1))
        for temp1 in sublists:
            for temp2 in sublists:
                trip_temp_tuple = [0] * n
                for j in range(len(temp1)):
                    trip_temp_tuple[temp1[j] - 1] = temp2[j]
                trip_temp = triple.BDTriple(trip_temp_tuple)
                if trip_temp.valid():
                    res.append(trip_temp)
    for trip1 in res:
        for trip2 in res:
            if trip1 != trip2 and are_iso(trip1, trip2):
                res.remove(trip2)
    return res

def are_iso(trip1: triple.BDTriple, trip2: triple.BDTriple):
    n = trip1.n
    g11 = trip1.g1
    g12 = trip1.g2
    g21 = trip2.g1
    g22 = trip2.g2

    if len(g11) != len(g21):
        return False

    if g22 == g11 and g21 == g12:
        return True

    diff = (g11[0] - g21[0]) % n
    for i in range(len(g11)):
        if (g11[i] - g21[i]) % n != diff or (g12[i] - g22[i]) % n != diff:
            return False
    return True