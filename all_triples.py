import random, triple
from itertools import permutations

def nonassociative_affine_triples(n : int) -> list[triple.BDTriple]:
    res = []
    temp = all_triples(n)
    print(len(temp))
    for trip in temp:
        temp = trip.tuple
        if set([x + 1 for x in range(n) if temp[x] != 0]).union(set(temp)) >= set([x + 1 for x in range(n)]) and trip.associative() is None:
            res.append(trip)
    return res

def is_affine(trip : triple.BDTriple) -> bool:
    n, temp = trip.n, trip.tuple
    return set([x + 1 for x in range(n) if temp[x] != 0]).union(set(temp)) >= set([x + 1 for x in range(n)])

def all_triples(n : int) -> list[triple.BDTriple]:
    res = []
    lst = [(x + 1) for x in range(n)]
    for i in range(n):
        sublists = []
        sublists.extend(permutations(lst, i))
        for temp1 in sublists:
            for temp2 in sublists:
                if set(temp1).union(set(temp2)) >= set([x + 1 for x in range(n)]) and temp1 != temp2:
                        trip_temp = triple.BDTriple(None, n=n, g1=temp1, g2=temp2)
                        if trip_temp.valid():
                            res.append(trip_temp)
    print(f"Number of all possible triples: {len(res)}")
    group = []
    i = 0
    while res:
        # print(f"The {i}th reduction leaves {len(res)} triples")
        trip = res[0]
        equiv_class = dihedral_action(trip)
        group.append(equiv_class)
        res = list(set(res) - set(equiv_class))
        i+=1

    mod_res = []
    for equiv_class in group:
        mod_res.append(equiv_class[0])
    return mod_res

# Generate the orbit of a triple under actions of the dihedral group
def dihedral_action(trip : triple.BDTriple) -> list[triple.BDTriple]:
    g1, g2, n = trip.g1, trip.g2, trip.n
    res = []

    def red(a, b):
        return (a - 1) % b + 1

    for i in range(n):
        g1_ri = [red(x + i, n) for x in g1]
        g2_ri = [red(x + i, n) for x in g2]
        res.append(triple.BDTriple(None, n=n, g1=g1_ri, g2=g2_ri))
        res.append(triple.BDTriple(None, n=n, g1=g2_ri, g2=g1_ri))
        if n % 2 == 0:
            g1_ris = [n + 1 - x for x in g1_ri]
            g2_ris = [n + 1 - x for x in g2_ri]
        else:
            g1_ris = [red(n + 2 - x, n) for x in g1_ri]
            g2_ris = [red(n + 2 - x, n) for x in g2_ri]
        res.append(triple.BDTriple(None, n=n, g1=g1_ris, g2=g2_ris))
        res.append(triple.BDTriple(None, n=n, g1=g2_ris, g2=g1_ris))

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