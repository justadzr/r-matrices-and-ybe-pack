import sympy as sp, triple, ybe, mat1, mat2, nonassoc_affine

# Need to declare the variables first
x = sp.Symbol("x")
h = sp.Symbol("h")
qn = sp.Symbol("qn")

triples = nonassoc_affine.nonassoc_affine_triples(7)
num = len(triples)
print(f"There are {num} nonassociative affine triples.")
 
incorrect = []
i = 0
for trip in triples:
    i += 1
    R = ybe.ggs_conjecture_rat_new(trip, x, qn)
    res = ybe.qybe_rat(R, x).simplify_rat()
    if str(res) != "0":
        print(res)
        incorrect.append(str(trip) + "\n")
    else:
        print(f"Triple {str(trip)}:\nchecked: {i}/{num}")

print("The incorrect triples are:")
print(incorrect)