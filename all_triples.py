import random, triple
from itertools import combinations

def nonassociative_triples(n: int) -> list[triple.BDTriple]:
    res = []
    temp = all_triples(n)
    for trip in temp:
        if trip.associative() is None:
            res.append(trip)
    return res

def all_triples(n : int) -> list[triple.BDTriple]:
    res = []
    lst = [(x + 1) for x in range(n)]
    for i in range(n):
        sublists = []
        sublists.extend(combinations(lst, i))
        for temp1 in sublists:
            for temp2 in sublists:
                trip_temp_tuple = [0] * n
                for j in range(len(temp1)):
                    trip_temp_tuple[temp1[j] - 1] = temp2[j]
                trip_temp = triple.BDTriple(trip_temp_tuple)
                if trip_temp.valid():
                    res.append(trip_temp)
    group = [[res[0]]]
    for trip in res:
        in_group = False
        for i in range(len(group)):
            in_group = False
            for trip_temp in group[i]:
                if trip != trip_temp and are_iso(trip, trip_temp):
                    in_group = True
                    group[i].append(trip)
                    break
            else:
                continue
            break
        if not in_group:
            group.append([trip])
    res = []
    for subgroup in group:
        res.append(subgroup[0])
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