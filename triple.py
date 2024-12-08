from operator import itemgetter
from itertools import groupby
from numpy import zeros
import sympy as sp
import mat2, mat1

# Write a BD triple of the affine untwisted sl(n) as a list of length n, indicating the image
# of \alpha_i under the transformation T: for instance
# [a_1, ..., a_n]
# means T(\alpha_1)=\alpha_{a_n}. Note if a_i is then \alpha_i is not in \Gamma_1.
# From now let 1, ..., n represents the roots of the affine untwisted sl(n). 
# Write 0 if the ith root is not in the domain of the transformation.

# The code here is not optimal at all: nested loops and unnecessary if/else are everywhere.
# But usually we don't need large triples (at least our matrix computation can't be too hard)
# since otherwise Sympy.simplify() would be too slow. So I just leave it as it is.

def left_end(n, connected):
    for i in range(len(connected)):
        if (connected[i] - 2 + n) % n + 1 not in connected:
            return connected[i]

def right_end(n, connected):
    for i in range(len(connected)):
        if connected[i] % n + 1 not in connected:
            return connected[i]

def dual_basis(basis):
    """
    Finds the dual basis of a given basis of simple roots of the weight space.

    Args:
        basis (list of int): A list of integers corresponding to a basis of α_i

    Returns:
        list of MatrixTensor1: The dual basis corresponding to the given basis of simple roots.
    """
    n = len(basis)
    # print(basis)
    N = n + 1
    aux = []
    for i in basis:
        temp = sp.MutableDenseNDimArray(zeros((N,)*2).astype(int))
        temp[(i-1) % N, (i-1) % N] = -sp.Integer(1)
        temp[i % N, i % N] = sp.Integer(1)
        aux.append(mat1.MatrixTensor1(N, temp))
    
    system = sp.zeros(n, n)
    for i, j in [(x, y) for x in range(n) for y in range(n)]:
        system[i, j] = 2 * int(i == j) - int(basis[j] == basis[i] % N + 1) \
            - int(basis[j] == (basis[i] -2 + N) % N + 1)
    inverse = system.inv()
    res = []
    for i in range(n):
        temp = sp.MutableDenseNDimArray(zeros((N,)*2).astype(int))
        temp_mat = mat1.MatrixTensor1(N, temp)
        for j in range(n):
            temp_mat += inverse[i, j] * aux[j]
        res.append(temp_mat)
    return res

class BDTriple:
    def __init__(self, triple):
        self.n = len(triple)
        self.g1 = [] + [(x+1) for x in range(self.n) if triple[x] != 0]
        self.g2 = [] + [triple[x-1] for x in self.g1]

    def valid(self) -> bool:
        """
        Checks if a Belavin-Drinfeld triple is valid.

        Returns:
            bool: True if the Belavin-Drinfeld triple is valid, False otherwise.
        """
        if len(set(self.g2)) != len(self.g1):
            print("Not bijective")
            return False

        n = self.n
        for i in range(len(self.g1)):
            if self.g1[i] % n + 1 in self.g1:
                ind = self.g1.index(self.g1[i] % n + 1)
                temp = abs(self.g2[i] - self.g2[ind])
                if temp != 1 and temp != n - 1:
                    print("Not orthogonal")
                    return False
            if (self.g1[i] - 2 + n) % n + 1 in self.g1:
                ind = self.g1.index((self.g1[i] - 2 + n) % n + 1)
                temp = abs(self.g2[i] - self.g2[ind])
                if temp != 1 and temp != n - 1:
                    print("Not orthogonal")
                    return False
            if self.g2[i] % n + 1 in self.g2:
                ind = self.g2.index(self.g2[i] % n + 1)
                temp = abs(self.g1[i] - self.g1[ind])
                if temp != 1 and temp != n - 1:
                    print("Not orthogonal")
                    return False
            if (self.g2[i] - 2 + n) % n + 1 in self.g2:
                ind = self.g2.index((self.g2[i] - 2 + n) % n + 1)
                temp = abs(self.g1[i] - self.g1[ind])
                if temp != 1 and temp != n - 1:
                    print("Not orthogonal")
                    return False

        for i in range(len(self.g1)):
            temp = self.g1[i]
            for k in range(self.n + 1):
                temp = self.T(temp)
                if temp == 0:
                    break
            if k >= self.n:
                print("Not nilpotent")
                return False
        return True
    
    def T(self, i):
        """
        Returns the index of T(α_i).

        Args:
            i (int): The index for which T(α_i) is to be retrieved.

        Returns:
            int: The index of T(α_i).
        """
        if i not in self.g1:
            return 0
        else:
            return self.g2[self.g1.index(i)]

    def __repr__(self):
        res = ""
        sub = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
        for i in range(len(self.g1)):
            res += f"\u03B1{self.g1[i]} -> \u03B1{self.g2[i]}\n".translate(sub)
        return res

    def connected_components_aux(self, G):
        n = self.n
        res = [None] * n
        m = len(G)
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
    
    def connected_components(self):
        """
        Returns all connected components of the Belavin-Drinfeld triple.

        Returns:
            list: A list of connected components.
        """
        return self.connected_components_aux(self.g1)
    
    def connected_components_img(self):
        """
        Returns all connected components in the image of the Belavin-Drinfeld triple.

        Returns:
            list: A list of connected components in the correct order.
        """
        components = self.connected_components()
        unclean = self.connected_components_aux(sorted(self.g2))
        temp = []
        if len(components) != len(unclean):
            print(unclean)
            raise Exception("Triple not valid: not orthogonal.")
        else: 
            for i in components:
                for j in unclean:
                    if self.T(i[0]) in j:
                        temp.append(j)
                        break
            return temp

    # 0 if preserves orientation, 1 o/w
    def orientation_check(self, connected, connected_img):
        n = self.n
        l1 = left_end(n, connected)
        l2 = left_end(n, connected_img)
        return int(l2 != self.T(l1))

    def associative(self) -> bool:
        """
        Checks if a Belavin-Drinfeld triple is associative by 
            checking the orientation of each connected component.

        Returns:
            bool: True if the Belavin-Drinfeld triple is associative, False otherwise.
        """
        if not self.valid():
            raise Exception("Triple not valid")
        else:
            components = self.connected_components()
            components_img = self.connected_components_img()
            n = self.n
            for i in range(len(components)):
                if 1 < len(components[i]) < n: 
                    # It is impossible for the length to be n. But I will leave it be.
                    if self.orientation_check(components[i], components_img[i]) == 1:
                        print(f"The orientation of" + \
                              " the connected component {components[i]} is reversed")
                        return False
            return True

    def choose_r0(self, only_return_s) -> mat2.MatrixTensor2:
        """
        Chooses one correct r_0 for a given Belavin-Drinfeld triple using the method described in 
        Schedler, T. (1999). "Verification of the GGS Conjecture for sl(n), n ≤ 12". 

        Args:
            only_return_s (bool): If True, the function returns s instead of r_0.

        Returns:
            MatrixTensor2: The chosen r_0 or s, depending on the value of `only_return_s`.

        References:
            Schedler, T. (1999). "Verification of the GGS Conjecture for sl(n), n ≤ 12". 
            Available at: https://arxiv.org/abs/9901079
        """
        if not self.valid():
            raise Exception("Triple not valid")
        else:
            n = self.n
            g1 = self.g1
            g2 = self.g2
            excl = -1
            for i in range(1, n+1):
                if i not in g1:
                    excl = i
                    break
            basis = [x for x in range(1, n+1) if x != excl]
            dual_basis_temp = dual_basis(basis)
            dual_basis_modified = []

            T_mat = sp.zeros(len(basis), len(basis))
            for i in range(len(basis)):
                if basis[i] in g1:
                    j = self.T(basis[i])
                    if j not in basis:
                        for p in range(len(basis)):
                            T_mat[p, i] = -1
                    else:
                        T_mat[basis.index(j), i] = 1
            transform = (sp.eye(len(basis)) - T_mat).inv()

            for i in range(len(basis)):
                temp = mat1.MatrixTensor1(n, sp.MutableDenseNDimArray(zeros((n,)*2).astype(int)))
                for j in range(len(basis)):
                    temp += transform[i, j] * dual_basis_temp[j]
                dual_basis_modified.append(temp)

            def aux_inner(triple, i, j):
                k = triple.T(i)
                l = triple.T(j)
                return 2 * (int(i == j) + int(k == j)) \
                    - 2 * int(l != 0) * (int(i == l) + int(k == l)) \
                    - int(i%n+1 == j) - int(i == j%n+1) - int(k%n+1 == j) - int(k == j%n+1) \
                    + int(l != 0) * (int(i%n+1 == l) + int(i == l%n+1) + \
                                     int(k%n+1 == l) + int(k == l%n+1))
            
            s = mat2.MatrixTensor2(n, sp.MutableDenseNDimArray(zeros((n,)*4).astype(int)), True)
            for i, j in [(x, y) for x in range(len(basis)) for y in range(len(basis))]:
                if basis[i] in g1:
                    s += sp.Rational(1, 2) * aux_inner(self, basis[i], basis[j]) * \
                        dual_basis_modified[i].tensor(dual_basis_modified[j])
                elif basis[j] in g1:
                    s -= sp.Rational(1, 2) * aux_inner(self, basis[j], basis[i]) * \
                        dual_basis_modified[i].tensor(dual_basis_modified[j])
            return s + int(not only_return_s) * sp.Rational(1, 2) * mat2.casimir(n)

    def C(self, alpha, beta):
        """
        Computes the orientation indicator C_{αβ} for a given Belavin-Drinfeld triple 
            using the method described in Schedler, T. (2000). "Proof of the GGS Conjecture". 

        Args:
            alpha (tuple): A tuple of two integers, representing (i, j).
            beta (tuple): A tuple of two integers, representing (k, l).

        Returns:
            int or None: The orientation indicator C_{αβ}, or None if undefined.

        References:
            Schedler, T. (2000). "Proof of the GGS Conjecture". 
            Available at: https://arxiv.org/abs/math/0009173
        """
        if not self.valid():
            raise Exception("Triple not valid")

        i, j = alpha
        k, l = beta
        n = self.n
        components_1 = self.connected_components()
        components_2 = self.connected_components_img()
        
        # First we check if \alpha is in \tilde{\Gamma_1} and if \beta is in \tilde{\Gamma_2}
        in_span_1 = False
        in_span_2 = False
        for connected_1 in components_1:
            if j in connected_1 and (i - 2 + n) % n + 1  in connected_1:
                in_span_1 = True
                break
        for connected_2 in components_2:
            if l in connected_2 and (k - 2 + n) % n + 1  in connected_2:
                in_span_2 = True
                break
        if not (in_span_1 and in_span_2):
            return None
        # The check here is redundant as we already have an orientation check
        # Yet this is probably better as we don't need to apply the check several times
        elif left_end(n, connected_2) == self.T(left_end(n, connected_1)):
            return 0
        else:
            return 1