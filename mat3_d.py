def hash(dim, i, j, k, l, p, q):
        return i * pow(dim, 5) + j * pow(dim, 4) + k * pow(dim, 3) + l * dim * dim + p * dim + q

class MatrixTensor3:
    def __init__(self, dim, coef):
        self.dim = dim
        self.coef = coef

    def __str__(self):
        dim = self.dim
        coef = self.coef
        result = ""
        for i, j in [(x, y) for x in range(dim) for y in range(dim)]:
            for k, l in [(x, y) for x in range(dim) for y in range(dim)]:
                for p, q in [(x, y) for x in range(dim) for y in range(dim)]:
                    temp = coef[hash(dim, i, j, k, l, p, q)]
                    sub = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
                    term = f"e{i+1}{j+1}⊗e{k+1}{l+1}⊗e{p+1}{q+1}".translate(sub)
                    if temp != 0:
                        result += "(" + str(temp) + ") * " + term + " + "
        return(result[:-3])

    def __add__(self, other):
        if self.dim == other.dim:
            dim = self.dim
            coef1 = self.coef
            coef2 = other.coef
            coef = coef1.copy()
            for i, j in [(x, y) for x in range(dim) for y in range(dim)]:
                for k, l in [(x, y) for x in range(dim) for y in range(dim)]:
                    for p, q in [(x, y) for x in range(dim) for y in range(dim)]:
                        coef[hash(dim, i, j, k, l, p, q)] = coef1[hash(dim, i, j, k, l, p, q)] + coef2[hash(dim, i, j, k, l, p, q)]
            return MatrixTensor3(dim, coef)
        else: 
            raise(ValueError, "Cannot add matrices of different dimensions")
    
    def __mul__(self, other):
        if self.dim == other.dim:
            dim = self.dim
            coef1 = self.coef
            coef2 = other.coef
            coef = coef1.copy()
            for i, j in [(x, y) for x in range(dim) for y in range(dim)]:
                for k, l in [(x, y) for x in range(dim) for y in range(dim)]:
                    for p, q in [(x, y) for x in range(dim) for y in range(dim)]:
                        coef[hash(dim, i, j, k, l, p, q)] = 0
                        for x, y, z in [(x, y, z) for x in range(dim) for y in range(dim) for z in range(dim)]:
                            coef[hash(dim, i, j, k, l, p, q)] += coef1[hash(dim, i, x, k, y, p, z)] * coef2[hash(dim, x, j, y, l, z, q)]
            return MatrixTensor3(dim, coef)
        else: 
            raise(ValueError, "Cannot add matrices of different dimensions")
