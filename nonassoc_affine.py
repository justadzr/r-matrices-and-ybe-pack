import triple, affine_diagram as ad
from itertools import chain, combinations
from anytree import Node
from time import time

def large_subsets(iterable : list):
    return chain.from_iterable(combinations(iterable, r) for 
                               r in range(-(-len(iterable)//2), len(iterable)))

def red(a, b):
    return (a - 1) % b + 1

def red_neg(a, b):
    return red(a, b) - b

def subdiagram(lst, n, left, right):
    if left > right:
        return lst[left-1:n] + lst[0:right]
    else:
        return lst[left-1:right]

def subdiagram_extend_one(lst, n, left, right):
    # THIS IS INCORRECT
    if left == 1:
        return [lst[-1]] + subdiagram(lst, n, left, right + 1)
    if right == n:
        return subdiagram(lst, n, left - 1, right) + [lst[0]]
    return subdiagram(lst, n, red(left - 1, n), red(right + 1, n))
            

def nonassoc_affine_triples(n : int) -> list[triple.BDTriple]:
    t = time()
    pretriples = []
    lists = large_subsets([x + 1 for x in range(n)])
    connected_components = []
    for lst in lists:
        connected_components += [ad.connected_components_aux(n, sorted(lst))]
  
    for g1_components in connected_components:
        max_depth = len(g1_components)
        g2_shift = [[] for _ in range(max_depth)] 
        roots = [[] for _ in range(max_depth)]
        where_g1_goes = [[] for _ in range(max_depth)]
        g1 = []
        g1_component_length = []
        g1_component_left_end = []

        max = 0
        max_ind = -1
        one_ind = -1
        for i in range(len(g1_components)):
            component = g1_components[i]
            g1_component_length += [len(component)]
            g1_component_left_end += [ad.left_end(n, component)]
            g1 += component
            if len(component) >= max:
                max = len(component)
                max_ind = i
            if g1_component_left_end[i] == 1:
                one_ind = i
        
        if one_ind > -1 and one_ind == max_ind:
            for shift in range(1, n):
                g2_shift[0] += [Node((shift, 0))]
                
                A = red(g1_component_left_end[0] + shift, n)
                B = red(g1_component_left_end[0] + g1_component_length[0] - 1 + shift, n)
                if A > B:
                    roots[0] += [Node([1] * B + [0] * (A - B - 1) + [1] * (n - A + 1))]
                else:
                    roots[0] += [Node([0] * (A - 1) + [1] * (B - A + 1) + [0] * (n - B))]

                if g1_component_length[0] > 1:
                    g2_shift[0] += [Node((shift, 1))]
                    if A > B:
                        roots[0] += [Node([1] * B + [0] * (A - B - 1) + [1] * (n - A + 1))]
                    else:
                        roots[0] += [Node([0] * (A - 1) + [1] * (B - A + 1) + [0] * (n - B))]

            depth = 0
            while depth < max_depth-1:
                depth += 1
                for i in range(len(g2_shift[depth-1])):
                    parent_node_shift = g2_shift[depth-1][i]
                    if parent_node_shift.name is not None:
                        parent_node_root = roots[depth-1][i]
                        for shift in range(1, n):
                            A = red(g1_component_left_end[depth] + shift, n)
                            B = red(g1_component_left_end[depth] + g1_component_length[depth] - 1 + shift, n)
                            if A > B:
                                current_diagram = subdiagram_extend_one(parent_node_root.name, n, A, B)
                                if 1 not in current_diagram:
                                    g2_shift[depth] += [Node((shift, 0), parent=parent_node_shift)]
                                    roots[depth] += [Node([1] * B + parent_node_root.name[B:A-1]
                                                        + [1] * (n - A + 1), parent=parent_node_root)]
                                    if g1_component_length[depth] > 1:
                                        g2_shift[depth] += [Node((shift, 1), parent=parent_node_shift)]
                                        roots[depth] += [Node([1] * B + parent_node_root.name[B:A-1]
                                                        + [1] * (n - A + 1), parent=parent_node_root)]
                                    
                                else:
                                    g2_shift[depth] += [Node(None, parent=parent_node_shift)]
                                    roots[depth] += [Node(None, parent=parent_node_root)]
                            else:
                                current_diagram = subdiagram_extend_one(parent_node_root.name, n, A, B)
                                if 1 not in current_diagram:
                                    g2_shift[depth] += [Node((shift, 0), parent=parent_node_shift)]
                                    roots[depth] += [Node(parent_node_root.name[0:A-1] + [1] * (B - A + 1)
                                                        + parent_node_root.name[B:n], parent=parent_node_root)]
                                    if g1_component_length[depth] > 1:
                                        g2_shift[depth] += [Node((shift, 1), parent=parent_node_shift)]
                                        roots[depth] += [Node(parent_node_root.name[0:A-1] + [1] * (B - A + 1)
                                                        + parent_node_root.name[B:n], parent=parent_node_root)]
                                    
                                else:
                                    g2_shift[depth] += [Node(None, parent=parent_node_shift)]
                                    roots[depth] += [Node(None, parent=parent_node_root)]
                    else:
                        g2_shift[depth] += [Node(None, parent=parent_node_shift)]
                        roots[depth] += [Node(None, parent=parent_node_root)]

            img_trees = g2_shift[-1]
            roots_trees = roots[-1]

            for x in range(len(img_trees)):
                tree = img_trees[x]
                roots_taken = roots_trees[x]
                g2 = []
                explicit = str(tree)[6:-2].split('/')[1:]
                explicit_roots = str(roots_taken)[6:-2].split('/')[-1][1:-1].split(', ')
                affine = True
                if explicit[-1] != "None":
                    for y in range(len(explicit_roots)):
                        if int(explicit_roots[y]) == 0 and (y + 1) not in g1:
                            affine = False
                            break
                    if affine:
                        for j in range(max_depth):
                            shift = int(explicit[j][1:-1].split(', ')[0])
                            orient = int(explicit[j][1:-1].split(', ')[1])
                            temp = list(map(lambda x : red(x + shift, n), g1_components[j]))
                            if orient == 0:
                                g2 += temp
                            elif orient == 1:
                                nonassoc = True
                                left = temp.index(ad.left_end(n, temp))
                                temp = (temp[left:] + temp[:left])[::-1]
                                g2 += temp[-left:] + temp[:-left]
                        pretrip = triple.BDTriple(None, n=n, g1=g1, g2=g2)
                        pretriples += [pretrip]
        else:
            continue

    t_new = time()
    print(f"Step 1 completed: all pretriples. Time passed: {t_new-t:.2f}s.")

    t = t_new
    group = []
    i = 0
    while pretriples:
        # print(f"The {i}th reduction leaves {len(pretriples)} triples")
        trip = pretriples[0]
        equiv_class = dihedral_action(trip)
        group.append(equiv_class)
        pretriples = [pretrip for pretrip in pretriples if pretrip not in equiv_class]
        i += 1
    t_new = time()
    print(f"Step 2 completed: equivalence classes. Time passed: {t_new-t:.2f}s")
    
    t = t_new
    res = []
    for equiv_class in group:
        pretrip = equiv_class[0]
        g1 = pretrip.g1
        g2 = pretrip.g2
        if pretrip.valid_ortho() and pretrip.affine_nonassociative():
            res.append(pretrip)
    t_new = time()
    print(f"Step 3 completed: nonassociative affine triple validation. Time passed: {t_new-t:.2f}s")

    return res

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


