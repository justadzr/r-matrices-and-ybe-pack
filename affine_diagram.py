from itertools import groupby
from operator import itemgetter

def left_end(n, connected):
    for i in range(len(connected)):
        if (connected[i] - 2 + n) % n + 1 not in connected:
            return connected[i]

def right_end(n, connected):
    for i in range(len(connected)):
        if connected[i] % n + 1 not in connected:
            return connected[i]
        
def connected_components_aux(n : int, G : list[int]):
        n = n
        res = [None] * n
        m = n
        key_1, key_n = None, None
        for k, g in groupby(enumerate(G), lambda x : (x[0] - x[1]) % m):
            temp = list(map(itemgetter(1), g))
            if 1 in temp:
                key_1 = k
            if n in temp:
                key_n = k
            if res[k] is None:
                res[k] = temp
            else:
                res[k] += temp
        if key_1 is not None and key_n is not None and key_1 != key_n:
            res[key_1] += res[key_n]
            res[key_n] = None
        return [x for x in res if x is not None]

def dihedral_equiv(n : int, lst1 : list[int], lst2 : list[int]):
    pass