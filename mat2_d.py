import mat3_d

def hash(dim, i, j, k, l):
        return i * pow(dim, 3) + j * dim * dim + k * dim + l

class MatrixTensor2:
    def __init__(self, dim, coef):
        self.dim = dim
        self.coef = coef

    def __str__(self):
        dim = self.dim
        coef = self.coef
        result = ""
        for i, j in [(x, y) for x in range(dim) for y in range(dim)]:
            for k, l in [(x, y) for x in range(dim) for y in range(dim)]:
                temp = coef[hash(dim, i, j, k, l)]
                sub = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
                term = f"e{i+1}{j+1}⊗e{k+1}{l+1}".translate(sub)
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
                    coef[hash(dim, i, j, k, l)] = coef1[hash(dim, i, j, k, l)] + coef2[hash(dim, i, j, k, l)]
            return MatrixTensor2(dim, coef)
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
                    coef[hash(dim, i, j, k, l)] = 0
                    for p, q in [(x, y) for x in range(dim) for y in range(dim)]:
                        coef[hash(dim, i, j, k, l)] += coef1[hash(dim, i, p, k, q)] * coef2[hash(dim, p, j, q, l)]
            return MatrixTensor2(dim, coef)
        else: 
            raise(ValueError, "Cannot add matrices of different dimensions")
    
    def swap(self):
        coef1 = self.coef
        dim = self.dim
        coef = coef1.copy()
        for i, j in [(x, y) for x in range(dim) for y in range(dim)]:
                for k, l in [(x, y) for x in range(dim) for y in range(dim)]:
                    coef[hash(dim, i, j, k, l)] = coef1[hash(dim, k, l, i, j)]
        return MatrixTensor2(dim, coef)

    def casimir(self):
        dim = self.dim
        coef = [0] * pow(dim, 4)
        for i, j in [(x, y) for x in range(dim) for y in range(dim)]:
            coef[hash(dim, i, j, j, i)] = 1
        return MatrixTensor2(dim, coef)
    
    def flip(self):
        return self.casimir() * self
    
    def to_matrixtensor3_12(self):
        coef1 = self.coef
        dim = self.dim
        coef = [0] * pow(dim, 6)
        for i in range(pow(dim, 4)):
            for j in range(dim):
                coef[dim*dim*i + dim*j + j] = coef1[i]
        return mat3_d.MatrixTensor3(dim, coef)

    def to_matrixtensor3_23(self):
        coef1 = self.coef
        dim = self.dim
        coef = [0] * pow(dim, 6)
        for i in range(pow(dim, 4)):
            for j in range(dim):
                coef[pow(dim, 4)*(dim+1)*j + i] = coef1[i]
        return mat3_d.MatrixTensor3(dim, coef)
        
    def to_matrixtensor3_13(self):
        coef1 = self.coef
        dim = self.dim
        coef = [0] * pow(dim, 6)
        for i in range(pow(dim, 4)):
            temp1 = i % (dim*dim)
            temp2 = i - temp1
            for j in range(dim):
                coef[dim*dim*temp2 + dim*dim*(dim+1)*j + temp1] = coef1[i]
        return mat3_d.MatrixTensor3(dim, coef)
