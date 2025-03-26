import sympy as sp, triple, ybe

# Need to declare the variables first
x = sp.Symbol("x")
h = sp.Symbol("h")
qn = sp.Symbol("qn")

# for n in range(4, 13):
#     with open(f"nonassociative-affine-triples-{n}.txt", "w") as f:
#         print([trip.tuple for trip in nonassoc_affine.nonassoc_affine_triples(n)], file=f)

for n in range(4, 13):
    with open(f"nonassociative-affine-triples-{n}.txt", "r") as f:
        unclean_triples = f.read()[2:-2].split('], [')
    triples = []
    for unclean in unclean_triples:
        t = str.maketrans({',': ' ', '[': ' ', ']': ' '})
        tuple = list(map(int, unclean.translate(t).split()))
        if tuple:
            triples += [triple.BDTriple(tuple)]

    with open(f"verification-record-{n}.txt", "w") as f1, open(f"incorrect.txt", "a+") as f2:
        incorrect = []
        i = 0
        num = len(triples)
        print(f"There are {len(triples)} nonassociative affine triples when n = {n}.", file=f1)
        if num > 0:
            for trip in triples:
                i += 1
                R = ybe.ggs_conjecture_rat(trip, x, qn)
                res = ybe.qybe_rat(R, x).simplify_rat()
                if str(res) != "0":
                    print(res)
                    incorrect.append(str(trip))
                else:
                    print(f"Triple {str(trip)}:\nchecked: {i}/{num}", file=f1)
            print(f"For n = {n}, our formula is incorrect for:", file=f2)
            print(incorrect, file=f2)