import sympy as sp, triple, ybe, multiprocessing as mp
import mat2

# Need to declare the variables first
x = sp.Symbol("x")
h = sp.Symbol("h")
qn = sp.Symbol("qn")
u = sp.Symbol("u")

# for n in range(4, 13):
#     with open(f"nonassociative-affine-triples-{n}.txt", "w") as f:
#         print([trip.tuple for trip in nonassoc_affine.nonassoc_affine_triples(n)], file=f)

def check_conjecture_for_triple(triples, x, qn):
    for trip in triples:
        R = ybe.ggs_conjecture_rat(trip, x, qn)
        res = ybe.qybe_rat(R, x)
        temp = res.simplify_rat()
        if str(temp) != "0":
            print(f"Triple {str(trip)}:\nfailed.")
            incorrect.append(str(trip))
        else:
            print(f"Triple {str(trip)}:\nchecked.")

def check_unitary_for_triple(triples, u, qn):
    for trip in triples:
        R = ybe.ggs_conjecture_rat(trip, sp.exp(sp.Rational(1, n) * u), qn)
        R_norm = (1/((qn ** sp.Rational(n, 2) - qn ** sp.Rational(-n, 2)) ** (-1) + (sp.exp(sp.Rational(1, 2) * u) - sp.exp(sp.Rational(-1, 2) * u)) ** (-1))) * R
        temp = (  R_norm * (  R_norm.subs(u, -u).swap()  )  ).simplify()
        if str(temp) != str(mat2.identity(trip.n)):
            print(f"Triple {str(trip)}:\nfailed.")
            incorrect.append(str(trip))
        else:
            print(f"Triple {str(trip)}:\nchecked.")

if __name__ == "__main__":
    for n in range(4, 9):
        with open(f"nonassociative-affine-triples-{n}.txt", "r") as f:
            unclean_triples = f.read()[2:-2].split('], [')
        triples = []
        for unclean in unclean_triples:
            t = str.maketrans({',': ' ', '[': ' ', ']': ' '})
            tuple = list(map(int, unclean.translate(t).split()))
            if tuple:
                triples += [triple.BDTriple(tuple)]

        with open(f"incorrect-unitary-{n}.txt", "w") as f:
            incorrect = []
            check_unitary_for_triple(triples, u, qn)
            # num = len(triples)
            # # print(f"There are {len(triples)} nonassociative affine triples when n = {n}.", file=f0)
            # if num > 0:
            #     k = (num if int(num / 8.0) == 0 else int(num / 8.0))
            #     jobs = []
            #     for i in range(0, num, k):
            #         seg_triples = triples[i:i+k]
            #         p = mp.Process(target=check_unitary_for_triple(triple, u, qn),
            #                         args=(seg_triples, u, qn))
            #         jobs += [p]
            #         p.start()

            #     for p in jobs:
            #         p.join()

            print(f"For n = {n}, our formula is incorrect for:", file=f)
            print(incorrect, file=f)