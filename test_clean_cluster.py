import sympy as sp, triple, ybe, multiprocessing as mp

# Need to declare the variables first
x = sp.Symbol("x")
h = sp.Symbol("h")
qn = sp.Symbol("qn")

# for n in range(4, 13):
#     with open(f"nonassociative-affine-triples-{n}.txt", "w") as f:
#         print([trip.tuple for trip in nonassoc_affine.nonassoc_affine_triples(n)], file=f)

def check_conjecture_for_triple(triples, x, qn, queue):
    incorrect = []
    for trip in triples:
        R = ybe.ggs_conjecture_rat(trip, x, qn)
        res = ybe.qybe_rat(R, x)
        temp = res.simplify_rat()
        if str(temp) != "0":
            print(f"Triple {str(trip)}:\nfailed.")
            incorrect.append(trip)
    queue.put(incorrect)

if __name__ == "__main__":
    for n in range(10, 13):
        with open(f"nonassociative-affine-triples-{n}.txt", "r") as f:
            unclean_triples = f.read()[2:-2].split('], [')
        triples = []
        for unclean in unclean_triples:
            t = str.maketrans({',': ' ', '[': ' ', ']': ' '})
            tuple = list(map(int, unclean.translate(t).split()))
            if tuple:
                triples += [triple.BDTriple(tuple)]
        
        incorrect = []
        num = len(triples)
        # print(f"There are {len(triples)} nonassociative affine triples when n = {n}.", file=f0)
        if num > 0:
            k = int(num / 128.0)
            processes = []
            incorrect_list = []
            queue = mp.Queue()

            for i in range(0, num, k):
                seg_triples = triples[i:i+k]
                p = mp.Process(target=check_conjecture_for_triple,
                                args=(seg_triples, x, qn, queue))
                processes += [p]
                p.start()
            for p in processes:
                incorrect = queue.get()
                incorrect_list += incorrect
            for p in processes:
                p.join()

            print(f"n={n} finished.")
            print(f"For n = {n}, our formula is incorrect for:\n{incorrect_list}.")
    print("All finished.")