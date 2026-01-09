import sympy as sp, mat1, mat2, mat3, triple
from numpy import zeros
from colorama import Fore, Style, init
init()

# Need to declare the variables first
(u1, u2, u3) = sp.symbols("u1:4")

def subs(num_or_symbol, x, u):
    if isinstance(num_or_symbol, sp.Expr):
        return num_or_symbol.subs(x, u)
    else:
        return num_or_symbol
    
def r_1213(r: mat2.MatrixTensor2, x: sp.Symbol):
    dim = r.dim
    coef = r.coef
    if r.ggs:
        coef_1213 = sp.MutableDenseNDimArray(zeros((dim,)*6).astype(int))
        for i, j in [(x, y) for x in range(dim) for y in range(dim)]:
            for k, l in [(x, y) for x in range(dim) for y in range(dim)]:
                for p in range(dim):
                    q = (i + k + p - j - l) % dim
                    coef_1213[i, j, k, l, p, q] = subs(coef[i, (k + i - l) % dim, k, l], x, u1 - u2) *\
                        subs(coef[(k + i - l) % dim, j, p, q], x, u1 - u3)
        return mat3.MatrixTensor3(dim, coef_1213, True)
    else:
        return r.to_matrixtensor3_12().subs(x, u1 - u2) * r.to_matrixtensor3_13().subs(x, u1 - u3)
    
def r_1312(r: mat2.MatrixTensor2, x: sp.Symbol):
    dim = r.dim
    coef = r.coef
    if r.ggs:
        coef_1312 = sp.MutableDenseNDimArray(zeros((dim,)*6).astype(int))
        for i, j in [(x, y) for x in range(dim) for y in range(dim)]:
            for k, l in [(x, y) for x in range(dim) for y in range(dim)]:
                for p in range(dim):
                    q = (i + k + p - j - l) % dim
                    coef_1312[i, j, k, l, p, q] = subs(coef[i, (l + j - k) % dim, p, q], x, u1 - u3) *\
                        subs(coef[(l + j - k) % dim, j, k, l], x, u1 - u2)
        return mat3.MatrixTensor3(dim, coef_1312, True)
    else:
        return r.to_matrixtensor3_13().subs(x, u1 - u3) * r.to_matrixtensor3_12().subs(x, u1 - u2)
    
def r_1223(r: mat2.MatrixTensor2, x: sp.Symbol):
    dim = r.dim
    coef = r.coef
    if r.ggs:
        coef_1223 = sp.MutableDenseNDimArray(zeros((dim,)*6).astype(int))
        for i, j in [(x, y) for x in range(dim) for y in range(dim)]:
            for k, l in [(x, y) for x in range(dim) for y in range(dim)]:
                for p in range(dim):
                    q = (i + k + p - j - l) % dim
                    coef_1223[i, j, k, l, p, q] = subs(coef[i, j, k, (k + i - j) % dim], x, u1 - u2) *\
                        subs(coef[(k + i - j) % dim, l, p, q], x, u2 - u3)
        return mat3.MatrixTensor3(dim, coef_1223, True)
    else:
        return r.to_matrixtensor3_12().subs(x, u1 - u2) * r.to_matrixtensor3_23().subs(x, u2 - u3)
    
def r_2312(r: mat2.MatrixTensor2, x: sp.Symbol):
    dim = r.dim
    coef = r.coef
    if r.ggs:
        coef_1213 = sp.MutableDenseNDimArray(zeros((dim,)*6).astype(int))
        for i, j in [(x, y) for x in range(dim) for y in range(dim)]:
            for k, l in [(x, y) for x in range(dim) for y in range(dim)]:
                for p in range(dim):
                    q = (i + k + p - j - l) % dim
                    coef_1213[i, j, k, l, p, q] = subs(coef[k, (l + j - i) % dim, p, q], x, u2 - u3) *\
                        subs(coef[i, j, (l + j - i) % dim, l], x, u1 - u2)
        return mat3.MatrixTensor3(dim, coef_1213, True)
    else:
        return r.to_matrixtensor3_23().subs(x, u2 - u3) * r.to_matrixtensor3_12().subs(x, u1 - u2)
    
def r_1323(r: mat2.MatrixTensor2, x: sp.Symbol):
    dim = r.dim
    coef = r.coef
    if r.ggs:
        coef_1323 = sp.MutableDenseNDimArray(zeros((dim,)*6).astype(int))
        for i, j in [(x, y) for x in range(dim) for y in range(dim)]:
            for k, l in [(x, y) for x in range(dim) for y in range(dim)]:
                for p in range(dim):
                    q = (i + k + p - j - l) % dim
                    coef_1323[i, j, k, l, p, q] = subs(coef[i, j, p, (i + p - j) % dim], x, u1 - u3) *\
                        subs(coef[k, l, (i + p - j) % dim, q], x, u2 - u3)
        return mat3.MatrixTensor3(dim, coef_1323, True)
    else:
        return r.to_matrixtensor3_13().subs(x, u1 - u3) * r.to_matrixtensor3_23().subs(x, u2 - u3)
    
def r_2313(r: mat2.MatrixTensor2, x: sp.Symbol):
    dim = r.dim
    coef = r.coef
    if r.ggs:
        coef_2313 = sp.MutableDenseNDimArray(zeros((dim,)*6).astype(int))
        for i, j in [(x, y) for x in range(dim) for y in range(dim)]:
            for k, l in [(x, y) for x in range(dim) for y in range(dim)]:
                for p in range(dim):
                    q = (i + k + p - j - l) % dim
                    coef_2313[i, j, k, l, p, q] = subs(coef[i, j, (k + p - l) % dim, q], x, u1 - u3) *\
                        subs(coef[k, l, p, (k + p - l) % dim], x, u2 - u3)
        return mat3.MatrixTensor3(dim, coef_2313, True)
    else:
        return r.to_matrixtensor3_23().subs(x, u2 - u3) * r.to_matrixtensor3_13().subs(x, u1 - u3)

# Note this function does not simplify the result.
def cybe(r: mat2.MatrixTensor2, x: sp.Symbol) -> mat3.MatrixTensor3:
    return r_1213(r, x) - r_1312(r, x) + r_1223(r, x) - r_2312(r, x) + r_1323(r, x) - r_2313(r, x)

def qybe1(R: mat2.MatrixTensor2, x: sp.Symbol) -> mat3.MatrixTensor3:
    if R.ggs:
        dim = R.dim
        c1 = R.coef
        c = sp.MutableDenseNDimArray(zeros((dim,)*6).astype(int))
        for i, j in [(x, y) for x in range(dim) for y in range(dim)]:
            for k, l in [(x, y) for x in range(dim) for y in range(dim)]:
                for m in range(dim):
                    n = (i + k + m - j - l) % dim
                    for p in range(dim):
                        q = (k + i + m - p - j) % dim
                        c[i, j, k, l, m, n] += subs(c1[i, (k + i - p) % dim, k, p], x, u1 - u2)\
                            * subs(c1[(k + i - p) % dim, j, m, q], x, u1 - u3)\
                            * subs(c1[p, l, q, n], x, u2 - u3)
        return mat3.MatrixTensor3(dim, c, True)
    else:
        R12 = R.subs(x, u1-u2).to_matrixtensor3_12()
        R13 = R.subs(x, u1-u3).to_matrixtensor3_13()
        R23 = R.subs(x, u2-u3).to_matrixtensor3_23()
        return R12 * R13 * R23
    
def qybe1_aux(R: mat2.MatrixTensor2, x: sp.Symbol, attention: list) -> mat3.MatrixTensor3:
    if R.ggs:
        dim = R.dim
        c1 = R.coef
        c = sp.MutableDenseNDimArray(zeros((dim,)*6).astype(int))
        for i, j in [(x, y) for x in range(dim) for y in range(dim)]:
            for k, l in [(x, y) for x in range(dim) for y in range(dim)]:
                for m in range(dim):
                    n = (i + k + m - j - l) % dim
                    print_index = attention == [i, j, k, l, m, n]
                    if print_index:
                            print("*****************1111********************")
                            print(f"The coefficient before {[i, j, k, l, m, n]} is given by the sum of")
                    for p in range(dim):
                        q = (k + i + m - p - j) % dim
                        A = subs(c1[i, (k + i - p) % dim, k, p], x, u1 - u2)
                        B = subs(c1[(k + i - p) % dim, j, m, q], x, u1 - u3)
                        C = subs(c1[p, l, q, n], x, u2 - u3)
                        D = A * B * C
                        if print_index and str(D) != "0":
                            print(f"A at {[i, (k + i - p) % dim, k, p]}")
                            print(A)
                            print(f"B at {[(k + i - p) % dim, j, m, q]}")
                            print(B)
                            print(f"C at {[p, l, q, n]}")
                            print(C)
                            print("With product:")
                            print(D)
                            print("=========================================")
                        c[i, j, k, l, m, n] += D
        return mat3.MatrixTensor3(dim, c, True)
    else:
        R12 = R.subs(x, u1-u2).to_matrixtensor3_12
        R13 = R.subs(x, u1-u3).to_matrixtensor3_13
        R23 = R.subs(x, u2-u3).to_matrixtensor3_23
        return R12 * R13 * R23

def qybe2(R: mat2.MatrixTensor2, x: sp.Symbol) -> mat3.MatrixTensor3:
    if R.ggs:
        dim = R.dim
        c1 = R.coef
        c = sp.MutableDenseNDimArray(zeros((dim,)*6).astype(int))
        for i, j in [(x, y) for x in range(dim) for y in range(dim)]:
            for k, l in [(x, y) for x in range(dim) for y in range(dim)]:
                for m in range(dim):
                    n = (i + k + m - j - l) % dim
                    for p in range(dim):
                        c[i, j, k, l, m, n] += subs(c1[k, p, m, (k + m - p) % dim], x, u2 - u3)\
                            * subs(c1[i, (j + l - p) % dim, (k + m - p) % dim, n], x, u1 - u3)\
                            * subs(c1[(j + l - p) % dim, j, p, l], x, u1 - u2)
        return mat3.MatrixTensor3(dim, c, True)
    else:
        R12 = R.subs(x, u1-u2).to_matrixtensor3_12()
        R13 = R.subs(x, u1-u3).to_matrixtensor3_13()
        R23 = R.subs(x, u2-u3).to_matrixtensor3_23()
        return R23 * R13 * R12

def qybe2_aux(R: mat2.MatrixTensor2, x: sp.Symbol, attention: list) -> mat3.MatrixTensor3:
    if R.ggs:
        dim = R.dim
        c1 = R.coef
        c = sp.MutableDenseNDimArray(zeros((dim,)*6).astype(int))
        for i, j in [(x, y) for x in range(dim) for y in range(dim)]:
            for k, l in [(x, y) for x in range(dim) for y in range(dim)]:
                for m in range(dim):
                    n = (i + k + m - j - l) % dim
                    print_index = attention == [i, j, k, l, m, n]
                    if print_index:
                            print("*****************2222********************")
                            print(f"The coefficient of {[i, j, k, l, m, n]} is given by the sum of")
                    for p in range(dim):
                        A = subs(c1[k, p, m, (k + m - p) % dim], x, u2 - u3)
                        B = subs(c1[i, (j + l - p) % dim, (k + m - p) % dim, n], x, u1 - u3)
                        C = subs(c1[(j + l - p) % dim, j, p, l], x, u1 - u2)
                        D = A * B * C
                        if print_index and str(D) != "0":
                            print(f"A at {[k, p, m, (k + m - p) % dim]}")
                            print(A)
                            print(f"B at {[i, (j + l - p) % dim, (k + m - p) % dim, n]}")
                            print(B)
                            print(f"C at {[(j + l - p) % dim, j, p, l]}")
                            print(C)
                            print("With product:")
                            print(D)
                            print("=========================================")
                        c[i, j, k, l, m, n] += D
        return mat3.MatrixTensor3(dim, c, True)
    else:
        R12 = R.subs(x, u1-u2).to_matrixtensor3_12
        R13 = R.subs(x, u1-u3).to_matrixtensor3_13
        R23 = R.subs(x, u2-u3).to_matrixtensor3_23
        return R23 * R13 * R12

# Note this function does not simplify the result.
def qybe(R: mat3.MatrixTensor3, x: sp.Symbol) -> mat3.MatrixTensor3:
    # print("Computing the QYBE of the given R-matrix......")
    return qybe1(R, x) - qybe2(R, x)

def check_continuous_datum(trip: triple.BDTriple, r0: mat2.MatrixTensor2):
    dim = r0.dim
    if r0 + r0.swap() - mat2.casimir(trip.n) != mat2.zero(dim):
        return False
    for i in trip.g1:
        if(r0.root_action_left(trip.T(i)) + r0.root_action_right(i) != mat1.zero(dim)):
            print(r0.root_action_left(trip.T(i)) + r0.root_action_right(i))
            return False
    return True

def to_trigonometric_solution(trip: triple.BDTriple, x: sp.Symbol, standard_part: bool) \
    -> mat2.MatrixTensor2:
    n = trip.n
    T = trip.T
    components = trip.connected_components()
    coef1 = sp.MutableDenseNDimArray(zeros((n,)*4).astype(int))
    coef2 = sp.MutableDenseNDimArray(zeros((n,)*4).astype(int))

    # First I choose r0
    r0 = trip.choose_r0(only_return_s = False)

    # Next I insert the standard part (minus its terms in the Cartan)
    for m in range(1, n):
        for i, j in [(x, y) for x in range(n) for y in range(n)]:
            if (i - j) % n == m:
                # The standard part
                coef1[i, j, j, i] += int(standard_part) * \
                    x ** m / (x ** n - 1)
                
    def in_one_component(i, j):
        for connected in components:
            p = (i - j) % n
            if set([red(j + k, n) for k in range(p)]).issubset(connected):
                return True
        return False
    
    def take_out_ind(symb):
        s = str(symb)
        if s[0] == "-":
            return int(s[-1]), int(s[2])
        else:
            return int(s[1]), int(s[-1])
    
    def red(a, b):
        return (a - 1) % b + 1

    # At last the nonstandard part
    e = sp.symbols(f"e1:{n + 1}")
    for m in range(1, n):
        for i, j in [(x, y) for x in range(n) for y in range(n)]:
                num = 0
                i_human = i + 1
                j_human = j + 1
                if (i - j) % n == m:
                    while in_one_component(i_human, j_human):
                        a = i_human - 1
                        b = j_human - 1
                        root = 0
                        
                        p = (a - b) % n
                        for q in range(p):
                            root += e[T(red(j_human + q, n)) % n] - e[T(red(j_human + q, n)) - 1]
                        num += 1
                        k_human, l_human = take_out_ind(root)
                        k, l = k_human - 1, l_human - 1
                        indicator = trip.C((i + 1, j + 1), (k_human, l_human), num)
                        if indicator is None:
                            break
                        else:
                            coef2[k, l, j, i] -= (-1) ** (indicator * (abs(a - b) - 1)) * x ** m
                            coef2[j, i, k, l] += (-1) ** (indicator * (abs(b - a) - 1)) * x ** (-m)
                            i_human, j_human = k_human, l_human
                    # print(f"i: {i_human}, j: {j_human}, k: {k_human}, l: {l_human}")
                    # print(f"indicator: {indicator}")
                    # print(f"Here counter is {counter}\n======================")

    return  mat2.MatrixTensor2(n, coef2, True) + int(standard_part) * \
        (mat2.MatrixTensor2(n, coef1, True) + (1 / (x ** n - 1)) * mat2.casimir_gl(n) + r0)

def ggs_conjecture_rat_passing_ord(trip: triple.BDTriple, x: sp.Symbol, q_nth: sp.Symbol) \
    -> mat2.MatrixTensor2:
    dic = {}
    n = trip.n
    T = trip.T
    coef1 = mat2.to_sparray(n, [0] * pow(n, 4))
    for i, j in [(x, y) for x in range(n) for y in range(n)]:
            if i != j:
                coef1[i, i, j, j] += 1 / (q_nth ** sp.Rational(n, 2) - q_nth ** sp.Rational(-n, 2))
    for i in range(n):
        coef1[i, i, i, i] += 1 / (q_nth ** n - 1) + 1 / (1 - x ** (-n))
    
    coef2 = mat2.to_sparray(n, [0] * pow(n, 4))
    for m_human in range(1, n):
        for i in range(n):
            coef2[(i + m_human) % n, i, i, (i + m_human) % n] += \
                x ** m_human / (x ** n - 1)
    std1, std2 = mat2.MatrixTensor2(n, coef1, True), mat2.MatrixTensor2(n, coef2, True)

    s = trip.choose_r0(only_return_s=True)
    # saving some computation power by not twisting `std2`
    standard_part = s.exp_rat(q_nth, n, True) * std1 * s.exp_rat(q_nth, n, True) + std2

    components = trip.connected_components()
    coef = mat2.to_sparray(n, [0] * pow(n, 4))
    coef_s = s.coef
    
    def in_one_component(i, j):
        for connected in components:
            p = (i - j) % n
            if set([red(j + k, n) for k in range(p)]).issubset(connected):
                return True
        return False

    def take_out_ind(symb):
        s = str(symb)
        if s[0] == "-":
            lst = s[1:].split(' + ')
            return int(lst[1][1:]), int(lst[0][1:])
        else:
            lst = s.split(' - ')
            return int(lst[0][1:]), int(lst[1][1:])

    def red(a, b):
        return (a - 1) % b + 1

    e = sp.symbols(f"e1:{n + 1}")
    for m in range(1, n):
        for i, j in [(x, y) for x in range(n) for y in range(n)]:
                num = 0
                i_human = i + 1
                j_human = j + 1
                if (i - j) % n == m:
                    record = []
                    while in_one_component(i_human, j_human):
                        a = i_human - 1
                        b = j_human - 1
                        root = 0

                        p = (a - b) % n
                        for q in range(p):
                            root += e[T(red(j_human + q, n)) % n] - e[T(red(j_human + q, n))-1]
                        k_human, l_human = take_out_ind(root)
                        num += 1
                        record.append((k_human, l_human))
                        k, l = k_human - 1, l_human - 1
                        indicator = trip.C((i + 1, j + 1), (k_human, l_human), num)
                        if indicator is None:
                            break
                        else:
                            # print(record)
                            root_length = (i - j) % n
                            passed = 0
                            half_passed = 0

                            root_left_to_beta = \
                                (red(k_human - root_length, n), red(l_human - root_length, n))
                            root_right_to_beta = \
                                (red(k_human + root_length, n), red(l_human + root_length, n))
                            root_left_to_alpha = \
                                (red(i + 1 - root_length, n), red(j + 1 - root_length, n))
                            root_right_to_alpha = \
                                (red(i + 1 + root_length, n), red(j + 1 + root_length, n))

                            if root_left_to_beta == (i + 1, j + 1):
                                half_passed += 1
                            if root_right_to_beta == (i + 1, j + 1):
                                half_passed += 1

                            # THIS IS INCORRECT
                            if root_left_to_beta in record:
                                ord_from_alpha_to_root = record.index(root_left_to_beta) + 1
                                ord_from_root_to_beta = num - ord_from_alpha_to_root
                                # print(ord_from_alpha_to_root)
                                if root_length > 1:
                                    # print((i + 1, j + 1), root_left_to_beta)
                                    # print(trip.C((i + 1, j + 1), root_left_to_beta, 
                                    #           ord_from_alpha_to_root))
                                    if trip.C(root_left_to_beta, (k_human, l_human), 
                                              ord_from_root_to_beta) == indicator:
                                        passed += 1
                                else:
                                    C = None
                                    if root_left_to_alpha in record:
                                        C = trip.C((i + 1, j), (k_human, l_human - 1), 
                                                   ord_from_alpha_to_root)
                                    if C is None and root_right_to_alpha in record:
                                        C = trip.C((i + 2, j + 1), (k_human, l_human - 1), 
                                                   ord_from_alpha_to_root)
                                    if C is not None and C == 0:
                                        passed += 1
                            
                            if root_right_to_beta in record:                                
                                ord_from_alpha_to_root = record.index(root_right_to_beta) + 1
                                ord_from_root_to_beta = num - ord_from_alpha_to_root
                                if root_length > 1:
                                    if trip.C(root_right_to_beta, (k_human, l_human), 
                                              ord_from_root_to_beta) == indicator:
                                        passed += 1
                                else:
                                    C = None
                                    if root_left_to_alpha in record:
                                        C = trip.C((i + 1, j), (k_human + 1, l_human), 
                                            ord_from_alpha_to_root)
                                    if C is None and root_right_to_alpha in record:
                                        C = trip.C((i + 2, j + 1), (k_human + 1, l_human), 
                                            ord_from_alpha_to_root)
                                    if C is not None and C == 0:
                                        passed += 1
                                                                
                            passing_order = sp.Rational(1, 2) * half_passed + passed
                            dic[((i+1, j+1), (k+1, l+1))] = passing_order
                            
                            ps = 1 - coef_s[i, i, k, k] - coef_s[j, j, l, l] + coef_s[i, i, l, l] + coef_s[j, j, k, k]
                            temp = sp.Rational(1, 2) * (passing_order 
                                                        - coef_s[i, i, l, l] - coef_s[j, j, k, k]
                                                        + indicator * (root_length - 1))
                            if True:
                                print(f"The passing order at alpha=({i+1},{j+1}) beta=({k+1}, {l+1}) for T^{num} is {passing_order}")
                                print("=============================================")
                            coef[k, l, j, i] -= (-1) ** (indicator * (root_length - 1)) * \
                                q_nth ** (n * temp) * x ** m
                            coef[j, i, k, l] += (-1) ** (indicator * (root_length - 1)) * \
                                q_nth ** (-n * temp) / (x ** m)
                            i_human, j_human = k_human, l_human
    return dic, (standard_part + mat2.MatrixTensor2(n, coef, True))

def ess_twist_paper(trip: triple.BDTriple, x: sp.Symbol, q_nth: sp.Symbol, focus_lst):
    print(f"Computing the twist using the inverse triple: {trip}")
    n = trip.n
    T = trip.T
    J_inv = [mat2.identity(n)] * n
    J21 = [mat2.identity(n)] * n
    J21_inv = [mat2.identity(n)] * n
    J = [mat2.identity(n)] * n
    components = trip.connected_components()
    
    def in_one_component(i, j):
        for connected in components:
            p = (i - j) % n
            if set([red(j + k, n) for k in range(p)]).issubset(connected):
                return True
        return False

    def take_out_ind(symb):
        s = str(symb)
        if s[0] == "-":
            lst = s[1:].split(' + ')
            return int(lst[1][1:]), int(lst[0][1:])
        else:
            lst = s.split(' - ')
            return int(lst[0][1:]), int(lst[1][1:])

    def red(a, b):
        return (a - 1) % b + 1
    
    qq = q_nth ** sp.Rational(n, 2)
    PL, PR = trip.passing_orders()
    e = sp.symbols(f"e1:{n + 1}")
    for m in range(1, n):
        for i, j in [(x, y) for x in range(n) for y in range(n)]:
                num = 0
                i_human = i + 1
                j_human = j + 1
                if (i - j) % n == m:
                    record = []
                    while in_one_component(i_human, j_human):
                        a = i_human - 1
                        b = j_human - 1
                        root = 0

                        p = (a - b) % n
                        for q in range(p):
                            root += e[T(red(j_human + q, n)) % n] - e[T(red(j_human + q, n))-1]
                        k_human, l_human = take_out_ind(root)
                        
                        num += 1
                        record.append((k_human, l_human))
                        k, l = k_human - 1, l_human - 1
                        indicator = trip.C((i + 1, j + 1), (k_human, l_human), num)
                        if indicator is None:
                            break
                        else:
                            a_passing_order = PR[((i+1, j+1), (k_human, l_human))] - PL[((i+1, j+1), (k_human, l_human))]
                            root_length = (i - j) % n
                            # print(f"The APS at alpha=({i+1},{j+1}) beta=({k+1}, {l+1}) for T^{num} is {a_passing_order} with orientation C = {indicator}")

                            coef_j21 = mat2.to_sparray(n, [0] * pow(n, 4))
                            coef_j21[i, j, l, k] = (-qq) ** (-indicator * (root_length - 1)) * qq ** (a_passing_order) * (qq - qq ** (-1)) * x ** (-m)
                            J21[num-1] = J21[num-1] + mat2.MatrixTensor2(n, coef_j21, True)

                            coef_j_inv = mat2.to_sparray(n, [0] * pow(n, 4))
                            coef_j_inv[l, k, i, j] = (-qq) ** (indicator * (root_length - 1)) * qq ** (a_passing_order) * (qq - qq ** (-1)) * x ** (m)
                            J_inv[num-1] = J_inv[num-1] - mat2.MatrixTensor2(n, coef_j_inv, True)

                            coef_j = mat2.to_sparray(n, [0] * pow(n, 4))
                            coef_j[l, k, i, j] = (-qq) ** (-indicator * (root_length - 1)) * qq ** (a_passing_order) * (qq - qq ** (-1)) * x ** (m)
                            J[num-1] = J[num-1] + mat2.MatrixTensor2(n, coef_j, True)

                            coef_j21 = mat2.to_sparray(n, [0] * pow(n, 4))
                            coef_j21[i, j, l, k] = (-qq) ** (indicator * (root_length - 1)) * qq ** (a_passing_order) * (qq - qq ** (-1)) * x ** (-m)
                            J21_inv[num-1] = J21_inv[num-1] - mat2.MatrixTensor2(n, coef_j21, True)


                            
                            # ps = 1 - coef_s[i, i, k, k] - coef_s[j, j, l, l] + coef_s[i, i, l, l] + coef_s[j, j, k, k]
                            # temp = sp.Rational(1, 2) * (passing_order 
                            #                             - coef_s[i, i, l, l] - coef_s[j, j, k, k]
                            #                             + indicator * (root_length - 1))
                            # if passing_order > 1:
                            #     print(f"For the triple: {trip.to_latex()}:")
                            #     print(f"The passing order at alpha=({i+1},{j+1}) beta=({k+1}, {l+1}) for T^{num} is {passing_order}")
                            #     print("=============================================")
                            # coef[k, l, j, i] -= (-1) ** (indicator * (root_length - 1)) * \
                            #     q_nth ** (n * temp) * x ** m
                            # coef[j, i, k, l] += (-1) ** (indicator * (root_length - 1)) * \
                            #     q_nth ** (-n * temp) / (x ** m)
                            i_human, j_human = k_human, l_human

    coef1 = mat2.to_sparray(n, [0] * pow(n, 4))
    for i, j in [(x, y) for x in range(n) for y in range(n)]:
            if i != j:
                coef1[i, i, j, j] += 1
    for i in range(n):
        coef1[i, i, i, i] += qq + x ** (n) * (qq - qq ** (-1)) / (1 - x ** (n))
    
    coef2 = mat2.to_sparray(n, [0] * pow(n, 4))
    for m_human in range(1, n):
        for i in range(n):
            coef2[i, (i + m_human) % n, (i + m_human) % n, i] += \
                 (qq - qq ** (-1)) * x ** (m_human) / (1 - x ** (n))
    std = mat2.MatrixTensor2(n, coef1, True) + mat2.MatrixTensor2(n, coef2, True)

    res = mat2.identity(n)
    for j in J21_inv[::-1]:
        res = res * j
    res = res * std
    for j in J:
        res = res * j
    
    def focus(lst, f):
        return lst[tuple(f)]
    
    res_schedler = mat2.identity(n)

    def attention(x, y, n, lst):
        for i, j in [(a, b) for a in range(n) for b in range(n)]:
            if x[lst[0], i, lst[2], j] != 0 and y[i, lst[1], j, lst[3]] not in [0, 1]:
                lst1 = [lst[0], i, lst[2], j]
                lst2 = [i, lst[1], j, lst[3]]
                temp1 = x[lst[0], i, lst[2], j]
                temp2 = y[i, lst[1], j, lst[3]]
                print(Fore.GREEN + f"At {list(map(lambda p: p+1, lst1))} and {list(map(lambda p: p+1, lst2))}, we multiply {temp1} and {temp2} to get {temp1 * temp2}" + Style.RESET_ALL)

    for j in J_inv[::-1]:
        print("*********************************")
        print(f"After factor {j-mat2.identity(n)} + I")
        attention(res_schedler.coef, j.coef, n, focus_lst)
        res_schedler = res_schedler * j
        if isinstance(focus(res_schedler.coef, focus_lst), sp.Expr):
            print(f"To be e{list(map(lambda p: p+1, focus_lst))}: {focus(res_schedler.coef, focus_lst).simplify()}")
        else:
            print(f"To be e{list(map(lambda p: p+1, focus_lst))}: {focus(res_schedler.coef, focus_lst)}")
    print("*********************************")
    # print(std.subs(x, 1/x).swap())
    print("*********************************")
    print("After standard")
    attention(res_schedler.coef, std.subs(x, 1/x).swap().coef, n, focus_lst)
    res_schedler = res_schedler * std.subs(x, 1/x).swap()
    if isinstance(focus(res_schedler.coef, focus_lst), sp.Expr):
        print(f"To be e{list(map(lambda p: p+1, focus_lst))}: {focus(res_schedler.coef, focus_lst).simplify()}")
    else:
            print(f"To be e{list(map(lambda p: p+1, focus_lst))}: {focus(res_schedler.coef, focus_lst)}")
    for j in J21:
        print("*********************************")
        print(f"After factor {j-mat2.identity(n)} + I")
        attention(res_schedler.coef, j.coef, n, focus_lst)
        res_schedler = res_schedler * j   
        if isinstance(focus(res_schedler.coef, focus_lst), sp.Expr):
            print(f"To be e{list(map(lambda p: p+1, focus_lst))}: {focus(res_schedler.coef, focus_lst).simplify()}")
        else:
            print(f"To be e{list(map(lambda p: p+1, focus_lst))}: {focus(res_schedler.coef, focus_lst)}")

    # print(res_schedler.simplify())

    # print("All J^-1 factors: ------------------------------------------------")
    # print([(j - mat2.identity(n)).simplify() for j in J_inv])
    # print("All J21 factors: ------------------------------------------------")
    # print([(j - mat2.identity(n)).simplify() for j in J21])
    # print("Inverse check: =================================================")
    # for i in range(n):
    #     print((J21_inv[i] * J21[i] - mat2.identity(n)).simplify())
    # print("J21 check: =================================================")
    # for i in range(n):
    #     print((J21[i] - J[i].subs(x, 1/x).swap()).simplify())
    # # print("Standard R")
    # # print(std)
    # print("Done printing in the method ess_twist ==========================")
                
    return J_inv, J21, res, res_schedler, std.subs(x, 1/x).swap()

def ess_twist(trip: triple.BDTriple, x: sp.Symbol, q_nth: sp.Symbol):
    print(f"Computing the twist using the inverse triple: {trip}")
    n = trip.n
    T = trip.T
    J_inv = [mat2.identity(n)] * n
    J21 = [mat2.identity(n)] * n
    J21_inv = [mat2.identity(n)] * n
    J = [mat2.identity(n)] * n
    components = trip.connected_components()
    
    def in_one_component(i, j):
        for connected in components:
            p = (i - j) % n
            if set([red(j + k, n) for k in range(p)]).issubset(connected):
                return True
        return False

    def take_out_ind(symb):
        s = str(symb)
        if s[0] == "-":
            lst = s[1:].split(' + ')
            return int(lst[1][1:]), int(lst[0][1:])
        else:
            lst = s.split(' - ')
            return int(lst[0][1:]), int(lst[1][1:])

    def red(a, b):
        return (a - 1) % b + 1
    
    qq = q_nth ** sp.Rational(n, 2)
    PL, PR = trip.passing_orders()
    e = sp.symbols(f"e1:{n + 1}")
    for m in range(1, n):
        for i, j in [(x, y) for x in range(n) for y in range(n)]:
                num = 0
                i_human = i + 1
                j_human = j + 1
                if (i - j) % n == m:
                    record = []
                    while in_one_component(i_human, j_human):
                        a = i_human - 1
                        b = j_human - 1
                        root = 0

                        p = (a - b) % n
                        for q in range(p):
                            root += e[T(red(j_human + q, n)) % n] - e[T(red(j_human + q, n))-1]
                        k_human, l_human = take_out_ind(root)
                        
                        num += 1
                        record.append((k_human, l_human))
                        k, l = k_human - 1, l_human - 1
                        indicator = trip.C((i + 1, j + 1), (k_human, l_human), num)
                        if indicator is None:
                            break
                        else:
                            a_passing_order = -PR[((i+1, j+1), (k_human, l_human))] + PL[((i+1, j+1), (k_human, l_human))]
                            root_length = (i - j) % n
                            # print(f"The APS at alpha=({i+1},{j+1}) beta=({k+1}, {l+1}) for T^{num} is {a_passing_order} with orientation C = {indicator}")

                            coef_j21 = mat2.to_sparray(n, [0] * pow(n, 4))
                            coef_j21[j, i, k, l] = (-qq) ** (-indicator * (root_length - 1)) * qq ** (a_passing_order) * (qq - qq ** (-1)) * x ** (-m)
                            J21[num-1] = J21[num-1] + mat2.MatrixTensor2(n, coef_j21, True)

                            coef_j_inv = mat2.to_sparray(n, [0] * pow(n, 4))
                            coef_j_inv[k, l, j, i] = (-qq) ** (indicator * (root_length - 1)) * qq ** (a_passing_order) * (qq - qq ** (-1)) * x ** (m)
                            J_inv[num-1] = J_inv[num-1] - mat2.MatrixTensor2(n, coef_j_inv, True)

                            coef_j = mat2.to_sparray(n, [0] * pow(n, 4))
                            coef_j[k, l, j, i] = (-qq) ** (-indicator * (root_length - 1)) * qq ** (a_passing_order) * (qq - qq ** (-1)) * x ** (m)
                            J[num-1] = J[num-1] + mat2.MatrixTensor2(n, coef_j, True)

                            coef_j21 = mat2.to_sparray(n, [0] * pow(n, 4))
                            coef_j21[j, i, k, l] = (-qq) ** (indicator * (root_length - 1)) * qq ** (a_passing_order) * (qq - qq ** (-1)) * x ** (-m)
                            J21_inv[num-1] = J21_inv[num-1] - mat2.MatrixTensor2(n, coef_j21, True)


                            
                            # ps = 1 - coef_s[i, i, k, k] - coef_s[j, j, l, l] + coef_s[i, i, l, l] + coef_s[j, j, k, k]
                            # temp = sp.Rational(1, 2) * (passing_order 
                            #                             - coef_s[i, i, l, l] - coef_s[j, j, k, k]
                            #                             + indicator * (root_length - 1))
                            # if passing_order > 1:
                            #     print(f"For the triple: {trip.to_latex()}:")
                            #     print(f"The passing order at alpha=({i+1},{j+1}) beta=({k+1}, {l+1}) for T^{num} is {passing_order}")
                            #     print("=============================================")
                            # coef[k, l, j, i] -= (-1) ** (indicator * (root_length - 1)) * \
                            #     q_nth ** (n * temp) * x ** m
                            # coef[j, i, k, l] += (-1) ** (indicator * (root_length - 1)) * \
                            #     q_nth ** (-n * temp) / (x ** m)
                            i_human, j_human = k_human, l_human

    coef1 = mat2.to_sparray(n, [0] * pow(n, 4))
    for i, j in [(x, y) for x in range(n) for y in range(n)]:
            if i != j:
                coef1[i, i, j, j] += 1
    for i in range(n):
        coef1[i, i, i, i] += qq + x ** (n) * (qq - qq ** (-1)) / (1 - x ** (n))
    
    coef2 = mat2.to_sparray(n, [0] * pow(n, 4))
    for m_human in range(1, n):
        for i in range(n):
            coef2[(i + m_human) % n, i, i, (i + m_human) % n] += \
                 (qq - qq ** (-1)) * x ** (m_human) / (1 - x ** (n))
    std = mat2.MatrixTensor2(n, coef1, True) + mat2.MatrixTensor2(n, coef2, True)

    res = mat2.identity(n)
    for j in J21_inv[::-1]:
        res = res * j
    res = res * std
    for j in J:
        res = res * j
    
    def focus(lst, f):
        return lst[*f]
    
    res_schedler = mat2.identity(n)

    def attention(x, y, n, lst):
        for i, j in [(a, b) for a in range(n) for b in range(n)]:
            if x[lst[0], i, lst[2], j] != 0 and y[i, lst[1], j, lst[3]] not in [0, 1]:
                lst1 = [lst[0], i, lst[2], j]
                lst2 = [i, lst[1], j, lst[3]]
                temp1 = x[lst[0], i, lst[2], j]
                temp2 = y[i, lst[1], j, lst[3]]
                print(f"At {list(map(lambda p: p+1, lst1))} and {list(map(lambda p: p+1, lst2))}, we multiply {temp1} and {temp2} to get {temp1 * temp2}")

    focus_lst = [1, 0, 2, 3]

    for j in J_inv[::-1]:        
        res_schedler = res_schedler * j
        print("*********************************")
        print(f"After factor {j-mat2.identity(n)} + I")
        print(f"To be e{list(map(lambda p: p+1, focus_lst))}: {(1 / (qq  ** (-1) - qq)*focus(res_schedler.coef, focus_lst)).simplify()}")
    res_schedler = res_schedler * std.subs(x, 1/x).swap()
    print("*********************************")
    print(std.subs(x, 1/x).swap())
    print("*********************************")
    print("After standard")
    print(f"To be e{list(map(lambda p: p+1, focus_lst))}: {(1 / (qq  ** (-1) - qq)*focus(res_schedler.coef, focus_lst)).simplify()}")
    for j in J21:
        res_schedler = res_schedler * j
        print("*********************************")
        print(f"After factor {j-mat2.identity(n)} + I")
        attention(res.coef, j.coef, n, focus_lst)
        print(f"To be e{list(map(lambda p: p+1, focus_lst))}: {(1 / (qq  ** (-1) - qq)*focus(res_schedler.coef, focus_lst)).simplify()}")

    print("All J^-1 factors: ------------------------------------------------")
    print([(j - mat2.identity(n)).simplify() for j in J_inv])
    print("All J21 factors: ------------------------------------------------")
    print([(j - mat2.identity(n)).simplify() for j in J21])
    # print("Inverse check: =================================================")
    # for i in range(n):
    #     print((J21_inv[i] * J21[i] - mat2.identity(n)).simplify())
    # print("J21 check: =================================================")
    # for i in range(n):
    #     print((J21[i] - J[i].subs(x, 1/x).swap()).simplify())
    # # print("Standard R")
    # # print(std)
    # print("Done printing in the method ess_twist ==========================")
                
    return J_inv, J21, res, res_schedler, std.subs(x, 1/x).swap()


def ess_twist_const(trip: triple.BDTriple, q_nth: sp.Symbol) \
    -> mat2.MatrixTensor2:
    n = trip.n
    T = trip.T
    
    J_inv = [mat2.identity(n)] * n
    J21 = [mat2.identity(n)] * n
    J21_inv = [mat2.identity(n)] * n
    J = [mat2.identity(n)] * n
    components = trip.connected_components()

    PL, PR = trip.passing_orders()
    
    def in_one_component(i, j):
        for connected in components:
            p = (i - j) % n
            if set([red(j + k, n) for k in range(p)]).issubset(connected):
                return True
        return False

    def take_out_ind(symb):
        s = str(symb)
        if s[0] == "-":
            lst = s[1:].split(' + ')
            return int(lst[1][1:]), int(lst[0][1:])
        else:
            lst = s.split(' - ')
            return int(lst[0][1:]), int(lst[1][1:])

    def red(a, b):
        return (a - 1) % b + 1
    
    qq = q_nth ** sp.Rational(n, 2)
    e = sp.symbols(f"e1:{n + 1}")
    for m in range(1, n):
        for i, j in [(x, y) for x in range(n) for y in range(n)]:
                num = 0
                i_human = i + 1
                j_human = j + 1
                if (i - j) % n == m:
                    record = []
                    while in_one_component(i_human, j_human):
                        a = i_human - 1
                        b = j_human - 1
                        root = 0

                        p = (a - b) % n
                        for q in range(p):
                            root += e[T(red(j_human + q, n)) % n] - e[T(red(j_human + q, n))-1]
                        k_human, l_human = take_out_ind(root)
                        
                        num += 1
                        record.append((k_human, l_human))
                        k, l = k_human - 1, l_human - 1
                        indicator = trip.C((i + 1, j + 1), (k_human, l_human), num)
                        if indicator is None:
                            break
                        else:
                            a_passing_order = -PR[((i+1, j+1), (k_human, l_human))] + PL[((i+1, j+1), (k_human, l_human))]
                            root_length = (i - j) % n

                            print(f"The APS at alpha=({i+1},{j+1}) beta=({k+1}, {l+1}) for T^{num} is {a_passing_order}")

                            coef_j21 = mat2.to_sparray(n, [0] * pow(n, 4))
                            coef_j21[j, i, k, l] = (-qq) ** (-indicator * (root_length - 1)) * qq ** (a_passing_order) * (qq - qq ** (-1)) 
                            J21[num-1] = J21[num-1] + mat2.MatrixTensor2(n, coef_j21, True)

                            coef_j_inv = mat2.to_sparray(n, [0] * pow(n, 4))
                            coef_j_inv[k, l, j, i] = (-qq) ** (indicator * (root_length - 1)) * qq ** (a_passing_order) * (qq - qq ** (-1))
                            J_inv[num-1] = J_inv[num-1] - mat2.MatrixTensor2(n, coef_j_inv, True)

                            coef_j = mat2.to_sparray(n, [0] * pow(n, 4))
                            coef_j[k, l, j, i] = (-qq) ** (-indicator * (root_length - 1)) * qq ** (a_passing_order) * (qq - qq ** (-1))
                            J[num-1] = J[num-1] + mat2.MatrixTensor2(n, coef_j, True)

                            coef_j21 = mat2.to_sparray(n, [0] * pow(n, 4))
                            coef_j21[j, i, k, l] = (-qq) ** (indicator * (root_length - 1)) * qq ** (a_passing_order) * (qq - qq ** (-1))
                            J21_inv[num-1] = J21_inv[num-1] - mat2.MatrixTensor2(n, coef_j21, True)


                            
                            # ps = 1 - coef_s[i, i, k, k] - coef_s[j, j, l, l] + coef_s[i, i, l, l] + coef_s[j, j, k, k]
                            # temp = sp.Rational(1, 2) * (passing_order 
                            #                             - coef_s[i, i, l, l] - coef_s[j, j, k, k]
                            #                             + indicator * (root_length - 1))
                            # if passing_order > 1:
                            #     print(f"For the triple: {trip.to_latex()}:")
                            #     print(f"The passing order at alpha=({i+1},{j+1}) beta=({k+1}, {l+1}) for T^{num} is {passing_order}")
                            #     print("=============================================")
                            # coef[k, l, j, i] -= (-1) ** (indicator * (root_length - 1)) * \
                            #     q_nth ** (n * temp) * x ** m
                            # coef[j, i, k, l] += (-1) ** (indicator * (root_length - 1)) * \
                            #     q_nth ** (-n * temp) / (x ** m)
                            i_human, j_human = k_human, l_human

    coef1 = mat2.to_sparray(n, [0] * pow(n, 4))
    for i, j in [(x, y) for x in range(n) for y in range(n)]:
        if i != j:
            coef1[i, i, j, j] += 1
    for i in range(n):
        coef1[i, i, i, i] += qq
    
    coef2 = mat2.to_sparray(n, [0] * pow(n, 4))
    for i, j in [(x, y) for x in range(n) for y in range(n)]:
        if i < j:
            coef2[i, j, j, i] = (qq - qq ** (-1))
    std = mat2.MatrixTensor2(n, coef1, True) + mat2.MatrixTensor2(n, coef2, True)

    res = mat2.identity(n)
    for j in J21_inv[::-1]:
        res = res * j
    res = res * std
    for j in J:
        res = res * j

    res_schedler = mat2.identity(n)
    for j in J_inv[::-1]:
        res_schedler = res_schedler * j
        
    res_schedler = res_schedler * std
    for j in J21:
        res_schedler = res_schedler * j
        
    
    # print("All J^-1 factors: ==============================================")
    # print(J_inv)
    # print("All J21 factors: =================================================")
    # print(J21)
    # print("Inverse check: =================================================")
    # for i in range(n):
    #     print((J_inv[i] * J[i] - mat2.identity(n)).simplify())
    # print("J21 check: =================================================")
    # for i in range(n):
    #     print((J21[i] - J[i].swap()).simplify())
    # # print("Standard R")
    # # print(std)
    # print("Done printing in the method ess_twist ==========================")
                
    return res, res_schedler


def ggs_conjecture_rat(trip: triple.BDTriple, x: sp.Symbol, q_nth: sp.Symbol) \
    -> mat2.MatrixTensor2:
    # print(f"Computing the formula using the original triple: {trip}")
    n = trip.n
    T = trip.T
    PL, PR = trip.passing_orders()
    coef1 = mat2.to_sparray(n, [0] * pow(n, 4))
    for i, j in [(x, y) for x in range(n) for y in range(n)]:
            if i != j:
                coef1[i, i, j, j] += 1 / (q_nth ** sp.Rational(n, 2) - q_nth ** sp.Rational(-n, 2))
    for i in range(n):
        coef1[i, i, i, i] += 1 / (q_nth ** n - 1) + 1 / (1 - x ** (-n))
    
    coef2 = mat2.to_sparray(n, [0] * pow(n, 4))
    for m_human in range(1, n):
        for i in range(n):
            coef2[(i + m_human) % n, i, i, (i + m_human) % n] += \
                x ** m_human / (x ** n - 1)
    std1, std2 = mat2.MatrixTensor2(n, coef1, True), mat2.MatrixTensor2(n, coef2, True)

    s = trip.choose_r0(only_return_s=True)
    # saving some computation power by not twisting `std2`
    standard_part = s.exp_rat(q_nth, n, True) * std1 * s.exp_rat(q_nth, n, True) + std2

    components = trip.connected_components()
    coef = mat2.to_sparray(n, [0] * pow(n, 4))
    coef_s = s.coef
    
    def in_one_component(i, j):
        for connected in components:
            p = (i - j) % n
            if set([red(j + k, n) for k in range(p)]).issubset(connected):
                return True
        return False

    def take_out_ind(symb):
        s = str(symb)
        if s[0] == "-":
            lst = s[1:].split(' + ')
            return int(lst[1][1:]), int(lst[0][1:])
        else:
            lst = s.split(' - ')
            return int(lst[0][1:]), int(lst[1][1:])

    def red(a, b):
        return (a - 1) % b + 1

    e = sp.symbols(f"e1:{n + 1}")
    for m in range(1, n):
        for i, j in [(x, y) for x in range(n) for y in range(n)]:
                num = 0
                i_human = i + 1
                j_human = j + 1
                if (i - j) % n == m:
                    record = []
                    while in_one_component(i_human, j_human):
                        a = i_human - 1
                        b = j_human - 1
                        root = 0

                        p = (a - b) % n
                        for q in range(p):
                            root += e[T(red(j_human + q, n)) % n] - e[T(red(j_human + q, n))-1]
                        k_human, l_human = take_out_ind(root)
                        
                        num += 1
                        record.append((k_human, l_human))
                        k, l = k_human - 1, l_human - 1
                        indicator = trip.C((i + 1, j + 1), (k_human, l_human), num)
                        if indicator is None:
                            break
                        else:
                            # print(record)
                            root_length = (i - j) % n
                            # passed = 0
                            # half_passed = 0

                            # root_left_to_beta = \
                            #     (red(k_human - root_length, n), red(l_human - root_length, n))
                            # root_right_to_beta = \
                            #     (red(k_human + root_length, n), red(l_human + root_length, n))
                            # root_left_to_alpha = \
                            #     (red(i + 1 - root_length, n), red(j + 1 - root_length, n))
                            # root_right_to_alpha = \
                            #     (red(i + 1 + root_length, n), red(j + 1 + root_length, n))

                            # if root_left_to_beta == (i + 1, j + 1):
                            #     half_passed += 1
                            # if root_right_to_beta == (i + 1, j + 1):
                            #     half_passed += 1

                            # # THIS IS INCORRECT
                            # if root_left_to_beta in record:
                            #     ord_from_alpha_to_root = record.index(root_left_to_beta) + 1
                            #     ord_from_root_to_beta = num - ord_from_alpha_to_root
                            #     # print(ord_from_alpha_to_root)
                            #     if root_length > 1:
                            #         # print((i + 1, j + 1), root_left_to_beta)
                            #         # print(trip.C((i + 1, j + 1), root_left_to_beta, 
                            #         #           ord_from_alpha_to_root))
                            #         if trip.C(root_left_to_beta, (k_human, l_human), 
                            #                   ord_from_root_to_beta) == indicator:
                            #             passed += 1
                            #     else:
                            #         C = None
                            #         if root_left_to_alpha in record:
                            #             C = trip.C((i + 1, j), (k_human, l_human - 1), 
                            #                        ord_from_alpha_to_root)
                            #         if C is None and root_right_to_alpha in record:
                            #             C = trip.C((i + 2, j + 1), (k_human, l_human - 1), 
                            #                        ord_from_alpha_to_root)
                            #         if C is not None and C == 0:
                            #             passed += 1
                            
                            # if root_right_to_beta in record:                                
                            #     ord_from_alpha_to_root = record.index(root_right_to_beta) + 1
                            #     ord_from_root_to_beta = num - ord_from_alpha_to_root
                            #     if root_length > 1:
                            #         if trip.C(root_right_to_beta, (k_human, l_human), 
                            #                   ord_from_root_to_beta) == indicator:
                            #             passed += 1
                            #     else:
                            #         C = None
                            #         if root_left_to_alpha in record:
                            #             C = trip.C((i + 1, j), (k_human + 1, l_human), 
                            #                 ord_from_alpha_to_root)
                            #         if C is None and root_right_to_alpha in record:
                            #             C = trip.C((i + 2, j + 1), (k_human + 1, l_human), 
                            #                 ord_from_alpha_to_root)
                            #         if C is not None and C == 0:
                            #             passed += 1
                                    
                            passing_order = PR[((i+1, j+1), (k_human, l_human))] + PL[((i+1, j+1), (k_human, l_human))]      
                            
                            # ps = 1 - coef_s[i, i, k, k] - coef_s[j, j, l, l] + coef_s[i, i, l, l] + coef_s[j, j, k, k]
                            temp = sp.Rational(1, 2) * (passing_order 
                                                        - coef_s[i, i, l, l] - coef_s[j, j, k, k]
                                                        + indicator * (root_length - 1))
                            
                            # print(f"The passing order at alpha=({i+1},{j+1}) beta=({k+1}, {l+1}) for T^{num} is {passing_order}")
                            coef[k, l, j, i] -= (-1) ** (indicator * (root_length - 1)) * \
                                q_nth ** (n * temp) * x ** m
                            coef[j, i, k, l] += (-1) ** (indicator * (root_length - 1)) * \
                                q_nth ** (-n * temp) / (x ** m)
                            i_human, j_human = k_human, l_human
    return (standard_part + mat2.MatrixTensor2(n, coef, True))

def ggs_conjecture_rat_new(trip: triple.BDTriple, x: sp.Symbol, q_nth: sp.Symbol) \
    -> mat2.MatrixTensor2:
    n = trip.n
    T = trip.T
    coef1 = mat2.to_sparray(n, [0] * pow(n, 4))
    for i, j in [(x, y) for x in range(n) for y in range(n)]:
            if i != j:
                coef1[i, i, j, j] += 1 / (q_nth ** sp.Rational(n, 2) - q_nth ** sp.Rational(-n, 2))
    for i in range(n):
        coef1[i, i, i, i] += 1 / (q_nth ** n - 1) + 1 / (1 - x ** (-n))
    
    coef2 = mat2.to_sparray(n, [0] * pow(n, 4))
    for m_human in range(1, n):
        for i in range(n):
            coef2[(i + m_human) % n, i, i, (i + m_human) % n] += \
                x ** m_human / (x ** n - 1)
    std1, std2 = mat2.MatrixTensor2(n, coef1, True), mat2.MatrixTensor2(n, coef2, True)

    s = trip.choose_r0(only_return_s=True)
    # saving some computation power by not twisting `std2`
    standard_part = s.exp_rat(q_nth, n, True) * std1 * s.exp_rat(q_nth, n, True) + std2

    components = trip.connected_components()
    coef = mat2.to_sparray(n, [0] * pow(n, 4))
    coef_s = s.coef
    
    def in_one_component(i, j):
        for connected in components:
            p = (i - j) % n
            if set([red(j + k, n) for k in range(p)]).issubset(connected):
                return True
        return False

    def take_out_ind(symb):
        s = str(symb)
        if s[0] == "-":
            lst = s[1:].split(' + ')
            return int(lst[1][1:]), int(lst[0][1:])
        else:
            lst = s.split(' - ')
            return int(lst[0][1:]), int(lst[1][1:])

    def red(a, b):
        return (a - 1) % b + 1

    e = sp.symbols(f"e1:{n + 1}")
    for m in range(1, n):
        for i, j in [(x, y) for x in range(n) for y in range(n)]:
                num = 0
                i_human = i + 1
                j_human = j + 1
                if (i - j) % n == m:
                    record = []
                    while in_one_component(i_human, j_human):
                        a = i_human - 1
                        b = j_human - 1
                        root = 0

                        p = (a - b) % n
                        for q in range(p):
                            root += e[T(red(j_human + q, n)) % n] - e[T(red(j_human + q, n))-1]
                        k_human, l_human = take_out_ind(root)
                        num += 1
                        record.append((k_human, l_human))
                        k, l = k_human - 1, l_human - 1
                        indicator = trip.C((i + 1, j + 1), (k_human, l_human), num)
                        if indicator is None:
                            break
                        else:
                            # print(record)
                            root_length = (i - j) % n
                            passed = 0
                            half_passed = 0

                            root_left_to_beta = \
                                (red(k_human - root_length, n), red(l_human - root_length, n))
                            root_right_to_beta = \
                                (red(k_human + root_length, n), red(l_human + root_length, n))
                            root_left_to_alpha = \
                                (red(i + 1 - root_length, n), red(j + 1 - root_length, n))
                            root_right_to_alpha = \
                                (red(i + 1 + root_length, n), red(j + 1 + root_length, n))

                            if root_left_to_beta == (i + 1, j + 1):
                                half_passed += 1
                            if root_right_to_beta == (i + 1, j + 1):
                                half_passed += 1

                            # THIS IS INCORRECT
                            if root_left_to_beta in record:
                                ord_from_alpha_to_root = record.index(root_left_to_beta) + 1
                                ord_from_root_to_beta = num - ord_from_alpha_to_root
                                # print(ord_from_alpha_to_root)
                                if root_length > 1:
                                    # print((i + 1, j + 1), root_left_to_beta)
                                    # print(trip.C((i + 1, j + 1), root_left_to_beta, 
                                    #           ord_from_alpha_to_root))
                                    if trip.C(root_left_to_beta, (k_human, l_human), 
                                              ord_from_root_to_beta) == indicator:
                                        passed += 1
                                else:
                                    C = None
                                    if root_left_to_alpha in record:
                                        C = trip.C((i + 1, j), (k_human, l_human - 1), 
                                                   ord_from_alpha_to_root)
                                    if C is None and root_right_to_alpha in record:
                                        C = trip.C((i + 2, j + 1), (k_human, l_human - 1), 
                                                   ord_from_alpha_to_root)
                                    if C is not None and C == 0:
                                        passed += 1
                            
                            if root_right_to_beta in record:                                
                                ord_from_alpha_to_root = record.index(root_right_to_beta) + 1
                                ord_from_root_to_beta = num - ord_from_alpha_to_root
                                if root_length > 1:
                                    if trip.C(root_right_to_beta, (k_human, l_human), 
                                              ord_from_root_to_beta) == indicator:
                                        passed += 1
                                else:
                                    C = None
                                    if root_left_to_alpha in record:
                                        C = trip.C((i + 1, j), (k_human + 1, l_human), 
                                            ord_from_alpha_to_root)
                                    if C is None and root_right_to_alpha in record:
                                        C = trip.C((i + 2, j + 1), (k_human + 1, l_human), 
                                            ord_from_alpha_to_root)
                                    if C is not None and C == 0:
                                        passed += 1
                                    
                            passing_order = sp.Rational(1, 2) * half_passed + passed        
                            
                            ps = 1 - coef_s[i, i, k, k] - coef_s[j, j, l, l] + coef_s[i, i, l, l] + coef_s[j, j, k, k]
                            temp = sp.Rational(1, 2) * (passing_order 
                                                        - coef_s[i, i, l, l] - coef_s[j, j, k, k]
                                                        + indicator * (root_length - 1))
                            if passing_order > 1:
                                print(f"For the triple: {trip.to_latex()}:")
                                print(f"The passing order at alpha=({i+1},{j+1}) beta=({k+1}, {l+1}) for T^{num} is {passing_order}")
                                print("=============================================")
                            coef[k, l, j, i] -= (-1) ** (indicator * (root_length - 1)) * \
                                q_nth ** (n * temp) * x ** m
                            coef[j, i, k, l] += (-1) ** (indicator * (root_length - 1)) * \
                                q_nth ** (-n * temp) / (x ** m)
                            i_human, j_human = k_human, l_human
    return (standard_part + mat2.MatrixTensor2(n, coef, True))


def ggs_conjecture_rat_t(trip: triple.BDTriple, x: sp.Symbol, q_nth: sp.Symbol) \
    -> mat2.MatrixTensor2:
    n = trip.n
    T = trip.T
    coef1 = mat2.to_sparray(n, [0] * pow(n, 4))
    for i, j in [(x, y) for x in range(n) for y in range(n)]:
            if i != j:
                coef1[i, i, j, j] += 1 / (q_nth ** sp.Rational(n, 2) - q_nth ** sp.Rational(-n, 2))
    for i in range(n):
        coef1[i, i, i, i] += 1 / (q_nth ** n - 1) + 1 / (1 - x ** (-n))
    
    coef2 = mat2.to_sparray(n, [0] * pow(n, 4))
    for m_human in range(1, n):
        for i in range(n):
            coef2[i, (i + m_human) % n, (i + m_human) % n, i] += \
                x ** m_human / (x ** n - 1)
    std1, std2 = mat2.MatrixTensor2(n, coef1, True), mat2.MatrixTensor2(n, coef2, True)

    s = trip.choose_r0(only_return_s=True)
    # saving some computation power by not twisting `std2`
    standard_part = s.exp_rat(q_nth, n, True) * std1 * s.exp_rat(q_nth, n, True) + std2

    components = trip.connected_components()
    coef = mat2.to_sparray(n, [0] * pow(n, 4))
    coef_s = s.coef
    
    def in_one_component(j, i):
        for connected in components:
            p = (j - i) % n
            if set([red(i + k, n) for k in range(p)]).issubset(connected):
                return True
        return False

    def take_out_ind(symb):
        s = str(symb)
        if s[0] == "-":
            lst = s[1:].split(' + ')
            return int(lst[1][1:]), int(lst[0][1:])
        else:
            lst = s.split(' - ')
            return int(lst[0][1:]), int(lst[1][1:])

    def red(a, b):
        return (a - 1) % b + 1
    
    def switch(tup):
        return (tup[1], tup[0])

    e = sp.symbols(f"e1:{n + 1}")
    for m in range(1, n):
        for i, j in [(x, y) for x in range(n) for y in range(n)]:
                num = 0
                i_human = i + 1
                j_human = j + 1
                if (j - i) % n == m:
                    record = []
                    while in_one_component(j_human, i_human):
                        a = i_human - 1
                        b = j_human - 1
                        root = 0

                        p = (b - a) % n
                        for q in range(p):
                            root += e[T(red(i_human + q, n))-1] - e[T(red(i_human + q, n)) % n]
                        k_human, l_human = take_out_ind(root)
                        num += 1
                        record.append((k_human, l_human))
                        k, l = k_human - 1, l_human - 1
                        indicator = trip.C((j + 1, i + 1), (l_human, k_human), num)
                        if indicator is None:
                            break
                        else:
                            # print(record)
                            root_length = (j - i) % n
                            passed = 0
                            half_passed = 0

                            root_left_to_beta = \
                                (red(k_human - root_length, n), red(l_human - root_length, n))
                            root_right_to_beta = \
                                (red(k_human + root_length, n), red(l_human + root_length, n))
                            root_left_to_alpha = \
                                (red(i + 1 - root_length, n), red(j + 1 - root_length, n))
                            root_right_to_alpha = \
                                (red(i + 1 + root_length, n), red(j + 1 + root_length, n))

                            if root_left_to_beta == (i + 1, j + 1):
                                half_passed += 1
                            if root_right_to_beta == (i + 1, j + 1):
                                half_passed += 1

                            # THIS IS INCORRECT
                            if root_left_to_beta in record:
                                ord_from_alpha_to_root = record.index(root_left_to_beta) + 1
                                ord_from_root_to_beta = num - ord_from_alpha_to_root
                                # print(ord_from_alpha_to_root)
                                if root_length > 1:
                                    # print((i + 1, j + 1), root_left_to_beta)
                                    # print(trip.C((i + 1, j + 1), root_left_to_beta, 
                                    #           ord_from_alpha_to_root))
                                    if trip.C(switch(root_left_to_beta), switch((k_human, l_human)), 
                                              ord_from_root_to_beta) == indicator:
                                        passed += 1
                                else:
                                    C = None
                                    if root_left_to_alpha in record:
                                        C = trip.C(switch((i, j+1)), switch((k_human-1, l_human)), 
                                                   ord_from_alpha_to_root)
                                    if C is None and root_right_to_alpha in record:
                                        C = trip.C(switch((i + 1, j + 2)), switch((k_human-1, l_human)), 
                                                   ord_from_alpha_to_root)
                                    if C is not None and C == 0:
                                        passed += 1
                            
                            if root_right_to_beta in record:                                
                                ord_from_alpha_to_root = record.index(root_right_to_beta) + 1
                                ord_from_root_to_beta = num - ord_from_alpha_to_root
                                if root_length > 1:
                                    if trip.C(switch(root_right_to_beta), switch((k_human, l_human)), 
                                              ord_from_root_to_beta) == indicator:
                                        passed += 1
                                else:
                                    C = None
                                    if root_left_to_alpha in record:
                                        C = trip.C(switch((i , j+ 1)), switch((k_human, l_human+1)), 
                                            ord_from_alpha_to_root)
                                    if C is None and root_right_to_alpha in record:
                                        C = trip.C(switch((i + 1, j + 2)), switch((k_human, l_human + 1)), 
                                            ord_from_alpha_to_root)
                                    if C is not None and C == 0:
                                        passed += 1
                                    
                            passing_order = sp.Rational(1, 2) * half_passed + passed        
                            
                            ps = 1 - coef_s[i, i, k, k] - coef_s[j, j, l, l] + coef_s[i, i, l, l] + coef_s[j, j, k, k]
                            temp = sp.Rational(1, 2) * (passing_order 
                                                        - coef_s[i, i, l, l] - coef_s[j, j, k, k]
                                                        + indicator * (root_length - 1))
                            if True:
                                print(i,j,k,l)
                                print(f"The passing order at alpha=({i+1},{j+1}) beta=({k+1}, {l+1}) is {passing_order}")
                                print(f"PS is P{ps}")
                                print(f"The s part is given by {coef_s[i, i, k, k] + coef_s[j, j, l, l]}")
                                print(f"The indicator is {indicator}")
                                print(f"The exponent is {temp}")
                                print("=============================================")
                            coef[k, l, j, i] -= (-1) ** (indicator * (root_length - 1)) * \
                                q_nth ** (n * temp) * x ** m
                            coef[j, i, k, l] += (-1) ** (indicator * (root_length - 1)) * \
                                q_nth ** (-n * temp) / (x ** m)
                            i_human, j_human = k_human, l_human
    return (standard_part + mat2.MatrixTensor2(n, coef, True))
    
def qybe1_rat(R: mat2.MatrixTensor2, x: sp.Symbol) -> mat3.MatrixTensor3:
    if R.ggs:
        dim = R.dim
        c1 = R.coef
        c = sp.MutableDenseNDimArray(zeros((dim,)*6).astype(int))
        for i, j in [(x, y) for x in range(dim) for y in range(dim)]:
            for k, l in [(x, y) for x in range(dim) for y in range(dim)]:
                for m in range(dim):
                    n = (i + k + m - j - l) % dim
                    for p in range(dim):
                        q = (k + i + m - p - j) % dim
                        c[i, j, k, l, m, n] += subs(c1[i, (k + i - p) % dim, k, p], x, u1 / u2)\
                            * subs(c1[(k + i - p) % dim, j, m, q], x, u1 / u3)\
                            * subs(c1[p, l, q, n], x, u2 / u3)
        return mat3.MatrixTensor3(dim, c, True)
    else:
        R12 = R.subs(x, u1 / u2).to_matrixtensor3_12()
        R13 = R.subs(x, u1 / u3).to_matrixtensor3_13()
        R23 = R.subs(x, u2 / u3).to_matrixtensor3_23()
        return R12 * R13 * R23
    
def qybe1_rat_no_eval(R: mat2.MatrixTensor2, x: sp.Symbol) -> mat3.MatrixTensor3:
    if R.ggs:
        dim = R.dim
        c1 = R.coef
        c = sp.MutableDenseNDimArray(zeros((dim,)*6).astype(int))
        for i, j in [(x, y) for x in range(dim) for y in range(dim)]:
            for k, l in [(x, y) for x in range(dim) for y in range(dim)]:
                for m in range(dim):
                    n = (i + k + m - j - l) % dim
                    for p in range(dim):
                        q = (k + i + m - p - j) % dim
                        c[i, j, k, l, m, n] += sp.Mul(sp.Mul(sp.Subs(c1[i, (k + i - p) % dim, k, p], x, u1 / u2, evaluate=False), 
                                                      sp.Subs(c1[(k + i - p) % dim, j, m, q], x, u1 / u3, evaluate=False), evaluate=False), 
                                                      sp.Subs(c1[p, l, q, n], x, u2 / u3, evaluate=False), evaluate=False)
        return mat3.MatrixTensor3(dim, c, True)
    else:
        R12 = R.subs(x, u1 / u2).to_matrixtensor3_12
        R13 = R.subs(x, u1 / u3).to_matrixtensor3_13
        R23 = R.subs(x, u2 / u3).to_matrixtensor3_23
        return R12 * R13 * R23
    
def qybe2_rat(R: mat2.MatrixTensor2, x: sp.Symbol) -> mat3.MatrixTensor3:
    if R.ggs:
        dim = R.dim
        c1 = R.coef
        c = sp.MutableDenseNDimArray(zeros((dim,)*6).astype(int))
        for i, j in [(x, y) for x in range(dim) for y in range(dim)]:
            for k, l in [(x, y) for x in range(dim) for y in range(dim)]:
                for m in range(dim):
                    n = (i + k + m - j - l) % dim
                    for p in range(dim):
                        c[i, j, k, l, m, n] += subs(c1[k, p, m, (k + m - p) % dim], x, u2 / u3)\
                            * subs(c1[i, (j + l - p) % dim, (k + m - p) % dim, n], x, u1 / u3)\
                            * subs(c1[(j + l - p) % dim, j, p, l], x, u1 / u2)
        return mat3.MatrixTensor3(dim, c, True)
    else:
        R12 = R.subs(x, u1 / u2).to_matrixtensor3_12()
        R13 = R.subs(x, u1 / u3).to_matrixtensor3_13()
        R23 = R.subs(x, u2 / u3).to_matrixtensor3_23()
        return R23 * R13 * R12
    
def qybe2_rat_no_eval(R: mat2.MatrixTensor2, x: sp.Symbol) -> mat3.MatrixTensor3:
    if R.ggs:
        dim = R.dim
        c1 = R.coef
        c = sp.MutableDenseNDimArray(zeros((dim,)*6).astype(int))
        for i, j in [(x, y) for x in range(dim) for y in range(dim)]:
            for k, l in [(x, y) for x in range(dim) for y in range(dim)]:
                for m in range(dim):
                    n = (i + k + m - j - l) % dim
                    for p in range(dim):
                        c[i, j, k, l, m, n] += sp.Mul(sp.Mul(sp.Subs(c1[k, p, m, (k + m - p) % dim], x, u2 / u3, evaluate=False), 
                                                             sp.Subs(c1[i, (j + l - p) % dim, (k + m - p) % dim, n], x, u1 / u3, evaluate=False), evaluate=False),
                                                             sp.Subs(c1[(j + l - p) % dim, j, p, l], x, u1 / u2, evaluate=False), evaluate=False)
        return mat3.MatrixTensor3(dim, c, True)
    else:
        R12 = R.subs(x, u1 / u2).to_matrixtensor3_12
        R13 = R.subs(x, u1 / u3).to_matrixtensor3_13
        R23 = R.subs(x, u2 / u3).to_matrixtensor3_23
        return R23 * R13 * R12

def qybe_rat(R: mat3.MatrixTensor3, x: sp.Symbol) -> mat3.MatrixTensor3:
    # print("Computing the QYBE of the given R-matrix......")
    return qybe1_rat(R, x) - qybe2_rat(R, x)

def qybe_rat_no_eval(R: mat3.MatrixTensor3, x: sp.Symbol) -> mat3.MatrixTensor3:
    # print("Computing the QYBE of the given R-matrix......")
    return qybe1_rat_no_eval(R, x) - qybe2_rat_no_eval(R, x)

def qybe1_rat_aux(R: mat2.MatrixTensor2, x: sp.Symbol, attention: list) -> mat3.MatrixTensor3:
    if R.ggs:
        dim = R.dim
        c1 = R.coef
        c = sp.MutableDenseNDimArray(zeros((dim,)*6).astype(int))
        for i, j in [(x, y) for x in range(dim) for y in range(dim)]:
            for k, l in [(x, y) for x in range(dim) for y in range(dim)]:
                for m in range(dim):
                    n = (i + k + m - j - l) % dim
                    print_index = attention == [i, j, k, l, m, n]
                    if print_index:
                            print("*****************1111********************")
                            print(f"The coefficient before {[i, j, k, l, m, n]} is given by the sum of")
                    for p in range(dim):
                        q = (k + i + m - p - j) % dim
                        A = subs(c1[i, (k + i - p) % dim, k, p], x, u1 / u2)
                        B = subs(c1[(k + i - p) % dim, j, m, q], x, u1 / u3)
                        C = subs(c1[p, l, q, n], x, u2 / u3)
                        D = A * B * C
                        if print_index and str(D) != "0":
                            print(f"A at {[i, (k + i - p) % dim, k, p]}")
                            print(A)
                            print(f"B at {[(k + i - p) % dim, j, m, q]}")
                            print(B)
                            print(f"C at {[p, l, q, n]}")
                            print(C)
                            print("With product:")
                            print(D)
                            print("=========================================")
                        c[i, j, k, l, m, n] += D
                    if print_index:
                        print(c[i, j, k, l, m, n].simplify())
        return mat3.MatrixTensor3(dim, c, True)
    else:
        R12 = R.subs(x, u1/u2).to_matrixtensor3_12()
        R13 = R.subs(x, u1/u3).to_matrixtensor3_13()
        R23 = R.subs(x, u2/u3).to_matrixtensor3_23()
        return R12 * R13 * R23

def qybe2_rat_aux(R: mat2.MatrixTensor2, x: sp.Symbol, attention: list) -> mat3.MatrixTensor3:
    if R.ggs:
        dim = R.dim
        c1 = R.coef
        c = sp.MutableDenseNDimArray(zeros((dim,)*6).astype(int))
        for i, j in [(x, y) for x in range(dim) for y in range(dim)]:
            for k, l in [(x, y) for x in range(dim) for y in range(dim)]:
                for m in range(dim):
                    n = (i + k + m - j - l) % dim
                    print_index = attention == [i, j, k, l, m, n]
                    if print_index:
                            print("*****************2222********************")
                            print(f"The coefficient of {[i, j, k, l, m, n]} is given by the sum of")
                    for p in range(dim):
                        A = subs(c1[k, p, m, (k + m - p) % dim], x, u2 / u3)
                        B = subs(c1[i, (j + l - p) % dim, (k + m - p) % dim, n], x, u1 / u3)
                        C = subs(c1[(j + l - p) % dim, j, p, l], x, u1 / u2)
                        D = A * B * C
                        if print_index and str(D) != "0":
                            print(f"A at {[k, p, m, (k + m - p) % dim]}")
                            print(A)
                            print(f"B at {[i, (j + l - p) % dim, (k + m - p) % dim, n]}")
                            print(B)
                            print(f"C at {[(j + l - p) % dim, j, p, l]}")
                            print(C)
                            print("With product:")
                            print(D)
                            print("=========================================")
                        c[i, j, k, l, m, n] += D
                    if print_index:
                        print(c[i, j, k, l, m, n].simplify())
        return mat3.MatrixTensor3(dim, c, True)
    else:
        R12 = R.subs(x, u1/u2).to_matrixtensor3_12()
        R13 = R.subs(x, u1/u3).to_matrixtensor3_13()
        R23 = R.subs(x, u2/u3).to_matrixtensor3_23()
        return R23 * R13 * R12

def qybe_rat_raw(R: mat3.MatrixTensor3, x: sp.Symbol) -> mat3.MatrixTensor3:
    R12 = R.subs(x, u1/u2).to_matrixtensor3_12()
    R13 = R.subs(x, u1/u3).to_matrixtensor3_13()
    R23 = R.subs(x, u2/u3).to_matrixtensor3_23()
    return R12 * R13 * R23 - R23 * R13 * R12

def to_constant_solution(trip: triple.BDTriple, standard_part: bool) \
    -> mat2.MatrixTensor2:
    n = trip.n
    T = trip.T
    components = trip.connected_components()
    coef1 = sp.MutableDenseNDimArray(zeros((n,)*4).astype(int))
    coef2 = sp.MutableDenseNDimArray(zeros((n,)*4).astype(int))

    # First I choose r0
    r0 = trip.choose_r0(only_return_s = False)

    # Next I insert the standard part (minus its terms in the Cartan)
    
    for i, j in [(x, y) for x in range(n) for y in range(n)]:
        if i < j:
            # The standard part
            coef1[i, j, j, i] += 1
        
                
    def in_one_component(i, j):
        for connected in components:
            p = (i - j) % n
            if set([red(j + k, n) for k in range(p)]).issubset(connected):
                return True
        return False
    
    def take_out_ind(symb):
        s = str(symb)
        if s[0] == "-":
            return int(s[-1]), int(s[2])
        else:
            return int(s[1]), int(s[-1])
    
    def red(a, b):
        return (a - 1) % b + 1

    # At last the nonstandard part
    e = sp.symbols(f"e1:{n + 1}")
    for m in range(1, n):
        for i, j in [(x, y) for x in range(n) for y in range(n)]:
                num = 0
                i_human = i + 1
                j_human = j + 1
                if (i - j) % n == m:
                    while in_one_component(i_human, j_human):
                        a = i_human - 1
                        b = j_human - 1
                        root = 0
                        
                        p = (a - b) % n
                        for q in range(p):
                            root += e[T(red(j_human + q, n)) % n] - e[T(red(j_human + q, n)) - 1]
                        num += 1
                        k_human, l_human = take_out_ind(root)
                        k, l = k_human - 1, l_human - 1
                        indicator = trip.C((i + 1, j + 1), (k_human, l_human), num)
                        if indicator is None:
                            break
                        else:
                            root_length = (a - b) % n
                            coef2[k, l, j, i] -= (-1) ** (indicator * (root_length - 1))
                            coef2[j, i, k, l] += (-1) ** (indicator * (root_length - 1))
                            i_human, j_human = k_human, l_human
                    # print(f"i: {i_human}, j: {j_human}, k: {k_human}, l: {l_human}")
                    # print(f"indicator: {indicator}")
                    # print(f"Here counter is {counter}\n======================")

    return mat2.MatrixTensor2(n, coef2, True) + int(standard_part) * (mat2.MatrixTensor2(n, coef1, True) + r0)

def ggs_conjecture_constant(trip: triple.BDTriple, q_nth: sp.Symbol) \
    -> mat2.MatrixTensor2:
    n = trip.n
    T = trip.T
    PL, PR = trip.passing_orders()
    qq = q_nth ** sp.Rational(n, 2)
    coef1 = mat2.to_sparray(n, [0] * pow(n, 4))
    for i, j in [(x, y) for x in range(n) for y in range(n)]:
        if i != j:
            coef1[i, i, j, j] += 1
    for i in range(n):
        coef1[i, i, i, i] += qq
    
    coef2 = mat2.to_sparray(n, [0] * pow(n, 4))
    for i, j in [(x, y) for x in range(n) for y in range(n)]:
        if i < j:
            coef2[i, j, j, i] = (qq - qq ** (-1))
    std1, std2 = mat2.MatrixTensor2(n, coef1, True), mat2.MatrixTensor2(n, coef2, True)

    s = trip.choose_r0(only_return_s=True)
    # saving some computation power by not twisting `std2`
    standard_part = s.exp_rat(q_nth, n, True) * std1 * s.exp_rat(q_nth, n, True) + std2

    components = trip.connected_components()
    coef = mat2.to_sparray(n, [0] * pow(n, 4))
    coef_s = s.coef
    
    def in_one_component(i, j):
        for connected in components:
            p = (i - j) % n
            if set([red(j + k, n) for k in range(p)]).issubset(connected):
                return True
        return False

    def take_out_ind(symb):
        s = str(symb)
        if s[0] == "-":
            lst = s[1:].split(' + ')
            return int(lst[1][1:]), int(lst[0][1:])
        else:
            lst = s.split(' - ')
            return int(lst[0][1:]), int(lst[1][1:])

    def red(a, b):
        return (a - 1) % b + 1

    e = sp.symbols(f"e1:{n + 1}")
    for m in range(1, n):
        for i, j in [(x, y) for x in range(n) for y in range(n)]:
                num = 0
                i_human = i + 1
                j_human = j + 1
                if (i - j) % n == m:
                    record = []
                    while in_one_component(i_human, j_human):
                        a = i_human - 1
                        b = j_human - 1
                        root = 0

                        p = (a - b) % n
                        for q in range(p):
                            root += e[T(red(j_human + q, n)) % n] - e[T(red(j_human + q, n))-1]
                        k_human, l_human = take_out_ind(root)
                        num += 1
                        record.append((k_human, l_human))
                        k, l = k_human - 1, l_human - 1
                        indicator = trip.C((i + 1, j + 1), (k_human, l_human), num)
                        if indicator is None:
                            break
                        else:
                            # print(record)
                            root_length = (i - j) % n
                                    
                            passing_order = PL[((i+1, j+1), (k_human, l_human))] + PR[((i+1, j+1), (k_human, l_human))]         

                            temp = sp.Rational(1, 2) * (passing_order 
                                                        - coef_s[i, i, l, l] - coef_s[j, j, k, k]
                                                        + indicator * (root_length - 1))
                            coef[k, l, j, i] -= (-1) ** (indicator * (root_length - 1)) * \
                                q_nth ** (n * temp) * (qq - qq ** (-1))
                            coef[j, i, k, l] += (-1) ** (indicator * (root_length - 1)) * \
                                q_nth ** (-n * temp) * (qq - qq ** (-1))
                            i_human, j_human = k_human, l_human
    return (standard_part + mat2.MatrixTensor2(n, coef, True))


def ggs_conjecture_constant_s(trip: triple.BDTriple, q_nth: sp.Symbol, s: mat2.MatrixTensor2) \
    -> mat2.MatrixTensor2:
    n = trip.n
    T = trip.T
    components = trip.connected_components()
    
    coef1 = mat2.to_sparray(n, [0] * pow(n, 4))
    for i in range(n):
        coef1[i, i, i, i] += 1
    std1 = mat2.MatrixTensor2(n, coef1, True)

    coef2 = mat2.to_sparray(n, [0] * pow(n, 4))
    for i, j in [(x, y) for x in range(n) for y in range(n)]:
        if i < j:
            # The standard part
            coef2[i, j, j, i] += 1
    std2 = mat2.MatrixTensor2(n, coef2, True)

    # saving some computation power by not twisting `std2`
    standard_part = 1 / (q_nth ** sp.Rational(n, 2)- q_nth ** sp.Rational(-n, 2)) * (
        (2 * s + std1).exp_rat(q_nth, n, q_exp_half=True)) + std2

    print(standard_part)

    components = trip.connected_components()
    coef = mat2.to_sparray(n, [0] * pow(n, 4))
    coef_s = s.coef        
                
    def in_one_component(i, j):
        for connected in components:
            p = (i - j) % n
            if set([red(j + k, n) for k in range(p)]).issubset(connected):
                return True
        return False
    
    def take_out_ind(symb):
        s = str(symb)
        if s[0] == "-":
            lst = s[1:].split(' + ')
            return int(lst[1][1:]), int(lst[0][1:])
        else:
            lst = s.split(' - ')
            return int(lst[0][1:]), int(lst[1][1:])
    
    def red(a, b):
        return (a - 1) % b + 1

    # At last the nonstandard part
    e = sp.symbols(f"e1:{n + 1}")
    for m in range(1, n):
        for i, j in [(x, y) for x in range(n) for y in range(n)]:
                num = 0
                i_human = i + 1
                j_human = j + 1
                if (i - j) % n == m:
                    while in_one_component(i_human, j_human):
                        a = i_human - 1
                        b = j_human - 1
                        root = 0
                        
                        p = (a - b) % n
                        for q in range(p):
                            root += e[T(red(j_human + q, n)) % n] - e[T(red(j_human + q, n)) - 1]
                        num += 1
                        k_human, l_human = take_out_ind(root)
                        k, l = k_human - 1, l_human - 1
                        indicator = trip.C((i + 1, j + 1), (k_human, l_human), num)
                        if indicator is None:
                            break
                        elif i > j and i_human > j_human:
                            ps = 1 - coef_s[i, i, k, k] - coef_s[j, j, l, l] + coef_s[i, i, l, l] + coef_s[j, j, k, k]
                            if False:
                                print(f"The passing order at alpha=({i+1},{j+1}) beta=({k+1}, {l+1}) is {ps}")
                                print(f"The s part is given by {coef_s[i, i, l, l] + coef_s[j, j, k, k]}")
                                print(f"The indicator is {indicator}")
                                print(f"The other s part is given by {coef_s[i, i, k, k] + coef_s[j, j, l, l]}")
                            root_length = a - b
                            temp = sp.Rational(1, 2) * (1 - coef_s[i, i, k, k] - coef_s[j, j, l, l] + indicator * (root_length - 1))
                            root_length = (a - b) % n
                            coef[k, l, j, i] -= (-1) ** (indicator * (root_length - 1)) *\
                                q_nth ** (n * temp)
                            coef[j, i, k, l] += (-1) ** (indicator * (root_length - 1))* \
                                q_nth ** (-n * temp)
                            i_human, j_human = k_human, l_human

    return mat2.MatrixTensor2(n, coef, True) + standard_part

def ggs_conjecture_constant_sch(trip: triple.BDTriple, q_nth: sp.Symbol) \
    -> mat2.MatrixTensor2:
    n = trip.n
    T = trip.T
    components = trip.connected_components()
    PL, PR = trip.passing_orders()
    
    coef1 = mat2.to_sparray(n, [0] * pow(n, 4))
    for i in range(n):
        coef1[i, i, i, i] += 1
    std1 = mat2.MatrixTensor2(n, coef1, True)

    coef2 = mat2.to_sparray(n, [0] * pow(n, 4))
    for i, j in [(x, y) for x in range(n) for y in range(n)]:
        if i < j:
            # The standard part
            coef2[i, j, j, i] += 1
    std2 = mat2.MatrixTensor2(n, coef2, True)

    s = trip.choose_r0(only_return_s=True)
    # saving some computation power by not twisting `std2`
    standard_part = 1 / (q_nth ** sp.Rational(n, 2)- q_nth ** sp.Rational(-n, 2)) * (
        (2 * s + std1).exp_rat(q_nth, n, q_exp_half=True)) + std2

    components = trip.connected_components()
    coef = mat2.to_sparray(n, [0] * pow(n, 4))
    coef_s = s.coef        
                
    def in_one_component(i, j):
        for connected in components:
            p = (i - j) % n
            if set([red(j + k, n) for k in range(p)]).issubset(connected):
                return True
        return False
    
    def take_out_ind(symb):
        s = str(symb)
        if s[0] == "-":
            lst = s[1:].split(' + ')
            return int(lst[1][1:]), int(lst[0][1:])
        else:
            lst = s.split(' - ')
            return int(lst[0][1:]), int(lst[1][1:])
    
    def red(a, b):
        return (a - 1) % b + 1

    # At last the nonstandard part
    e = sp.symbols(f"e1:{n + 1}")
    for m in range(1, n):
        for i, j in [(x, y) for x in range(n) for y in range(n)]:
                num = 0
                i_human = i + 1
                j_human = j + 1
                if (i - j) % n == m:
                    while in_one_component(i_human, j_human):
                        a = i_human - 1
                        b = j_human - 1
                        root = 0
                        
                        p = (a - b) % n
                        for q in range(p):
                            root += e[T(red(j_human + q, n)) % n] - e[T(red(j_human + q, n)) - 1]
                        num += 1
                        k_human, l_human = take_out_ind(root)
                        k, l = k_human - 1, l_human - 1
                        indicator = trip.C((i + 1, j + 1), (k_human, l_human), num)
                        if indicator is None:
                            break
                        elif i > j and i_human > j_human:
                            ps = 1 - coef_s[i, i, k, k] - coef_s[j, j, l, l] + coef_s[i, i, l, l] + coef_s[j, j, k, k]
                            if False:
                                print(f"The passing order at alpha=({i+1},{j+1}) beta=({k+1}, {l+1}) is {ps}")
                                print(f"The s part is given by {coef_s[i, i, l, l] + coef_s[j, j, k, k]}")
                                print(f"The indicator is {indicator}")
                                print(f"The other s part is given by {coef_s[i, i, k, k] + coef_s[j, j, l, l]}")
                            root_length = a - b
                            temp = sp.Rational(1, 2) * (1 - coef_s[i, i, k, k] - coef_s[j, j, l, l] + indicator * (root_length - 1))
                            root_length = (a - b) % n
                            coef[k, l, j, i] -= (-1) ** (indicator * (root_length - 1)) *\
                                q_nth ** (n * temp)
                            coef[j, i, k, l] += (-1) ** (indicator * (root_length - 1))* \
                                q_nth ** (-n * temp)
                            i_human, j_human = k_human, l_human

    return mat2.MatrixTensor2(n, coef, True) + standard_part
    
def ggs_conjecture_constant_sch_s(trip: triple.BDTriple, q_nth: sp.Symbol, s: mat2.MatrixTensor2) \
    -> mat2.MatrixTensor2:
    n = trip.n
    T = trip.T
    components = trip.connected_components()
    
    coef1 = mat2.to_sparray(n, [0] * pow(n, 4))
    for i in range(n):
        coef1[i, i, i, i] += 1
    std1 = mat2.MatrixTensor2(n, coef1, True)

    coef2 = mat2.to_sparray(n, [0] * pow(n, 4))
    for i, j in [(x, y) for x in range(n) for y in range(n)]:
        if i < j:
            # The standard part
            coef2[i, j, j, i] += 1
    std2 = mat2.MatrixTensor2(n, coef2, True)

    # saving some computation power by not twisting `std2`
    standard_part = 1 / (q_nth ** sp.Rational(n, 2)- q_nth ** sp.Rational(-n, 2)) * (
        (2 * s + std1).exp_rat(q_nth, n, q_exp_half=True)) + std2

    components = trip.connected_components()
    coef = mat2.to_sparray(n, [0] * pow(n, 4))
    coef_s = s.coef        
                
    def in_one_component(i, j):
        for connected in components:
            p = (i - j) % n
            if set([red(j + k, n) for k in range(p)]).issubset(connected):
                return True
        return False
    
    def take_out_ind(symb):
        s = str(symb)
        if s[0] == "-":
            lst = s[1:].split(' + ')
            return int(lst[1][1:]), int(lst[0][1:])
        else:
            lst = s.split(' - ')
            return int(lst[0][1:]), int(lst[1][1:])
    
    def red(a, b):
        return (a - 1) % b + 1

    # At last the nonstandard part
    e = sp.symbols(f"e1:{n + 1}")
    for m in range(1, n):
        for i, j in [(x, y) for x in range(n) for y in range(n)]:
                num = 0
                i_human = i + 1
                j_human = j + 1
                if (i - j) % n == m:
                    while in_one_component(i_human, j_human):
                        a = i_human - 1
                        b = j_human - 1
                        root = 0
                        
                        p = (a - b) % n
                        for q in range(p):
                            root += e[T(red(j_human + q, n)) % n] - e[T(red(j_human + q, n)) - 1]
                        num += 1
                        k_human, l_human = take_out_ind(root)
                        k, l = k_human - 1, l_human - 1
                        indicator = trip.C((i + 1, j + 1), (k_human, l_human), num)
                        if indicator is None:
                            break
                        elif i > j and i_human > j_human:
                            ps = 1 - coef_s[i, i, k, k] - coef_s[j, j, l, l] + coef_s[i, i, l, l] + coef_s[j, j, k, k]
                            if False:
                                print(f"The passing order at alpha=({i+1},{j+1}) beta=({k+1}, {l+1}) is {ps}")
                                print(f"The s part is given by {coef_s[i, i, l, l] + coef_s[j, j, k, k]}")
                                print(f"The indicator is {indicator}")
                                print(f"The other s part is given by {coef_s[i, i, k, k] + coef_s[j, j, l, l]}")
                            root_length = a - b
                            temp = sp.Rational(1, 2) * (1 - coef_s[i, i, k, k] - coef_s[j, j, l, l] + indicator * (root_length - 1))
                            root_length = (a - b) % n
                            coef[k, l, j, i] -= (-1) ** (indicator * (root_length - 1)) *\
                                q_nth ** (n * temp)
                            coef[j, i, k, l] += (-1) ** (indicator * (root_length - 1))* \
                                q_nth ** (-n * temp)
                            i_human, j_human = k_human, l_human

    return mat2.MatrixTensor2(n, coef, True) + standard_part