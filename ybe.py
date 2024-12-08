import mat2
import triple
import sympy as sp
from numpy import zeros

# Need to declare the variables first
(u1, u2, u3) = sp.symbols("u1:4")

# I should write these in the class MatrixTensor2, 
# but in that case we can't substitute explicitly.
def cybe(r, x):
    r12 = r.subs(x, u1-u2).to_matrixtensor3_12()
    r13 = r.subs(x, u1-u3).to_matrixtensor3_13()
    r23 = r.subs(x, u2-u3).to_matrixtensor3_23()
    return r12.comm(r13) + r12.comm(r23) + r13.comm(r23)

def cybe_eval(r, x, u1_val, u2_val, u3_val, x_val):
    r12 = r.subs(x, u1-u2).to_matrixtensor3_12()
    r13 = r.subs(x, u1-u3).to_matrixtensor3_13()
    r23 = r.subs(x, u2-u3).to_matrixtensor3_23()
    return r12.comm(r13) + r12.comm(r23) + r13.comm(r23)

def qybe(R, x):
    R12 = R.subs(x, u1-u2).to_matrixtensor3_12
    R13 = R.subs(x, u1-u3).to_matrixtensor3_13
    R23 = R.subs(x, u2-u3).to_matrixtensor3_23
    return R12 * R13 * R23 - R23 * R13 * R12

def check_continuous_datum(trip, r0: mat2.MatrixTensor2):
    if str(r0 + r0.swap() - mat2.casimir(trip.n)) != "0":
        return False
    for i in trip.g1:
        if(str(r0.root_action_left(trip.T(i)) + r0.root_action_right(i)) != "0"):
            print(r0.root_action_left(trip.T(i)) + r0.root_action_right(i))
            return False
    return True

def to_trigonometric_solution(triple: triple.BDTriple, x: sp.Symbol, standard_part: bool) \
    -> mat2.MatrixTensor2:
    n = triple.n
    T = triple.T
    components = triple.connected_components()
    components_img = triple.connected_components_img()
    coef1 = sp.MutableDenseNDimArray(zeros((n,)*4).astype(int))
    coef2 = sp.MutableDenseNDimArray(zeros((n,)*4).astype(int))

    # First I choose r0
    r0 = triple.choose_r0(only_return_s = False)

    # Next I insert the standard part (minus its terms in the Cartan)
    for m in range(1, n):
        for i, j in [(x, y) for x in range(n) for y in range(n)]:
            if (i - j) % n == m:
                # The standard part
                coef1[i, j, j, i] += int(standard_part) * \
                    sp.exp(sp.Rational(m, n) * x) / (sp.exp(x) - 1)

    # # At last the nonstandard part
    # for connected in components:
    #     for m in range(1, n):
    #         for i, j in [(x, y) for x in connected for y in connected]:
    #             if (i - j) % n == m

    for m in range(1, n):
        for i, j in [(x, y) for x in range(n) for y in range(n)]:
                i_human = i + 1
                j_human = j + 1
                # The nonstandard part
                # First apply T (if applicable)
                if i_human > j_human:
                    k_human, l_human = (T(i_human - 1) % n + 1, T(j_human))
                    k = k_human - 1
                    l = l_human - 1
                    indicator = triple.C((i_human, j_human), (k_human, l_human))
                else:
                    l_human, k_human = (T(i_human - 1) % n + 1, T(j_human))
                    k = k_human - 1
                    l = l_human - 1
                    indicator = triple.C((i_human, j_human), (k_human, l_human))
                # print(f"m = {m}")
                # print(f"i: {i_human}, j: {j_human}, k: {k_human}, l: {l_human}")
                # print(f"indicator: {indicator}")
                # counter = 1
                # print(f"Here counter is {counter}\n======================")

                # When the previous T is applicable
                while indicator is not None:
                    # Add one nonstandard term
                    coef2[k, l, j, i] -= (-1) ** (indicator * (abs(i - j) - 1)) * \
                        sp.exp(sp.Rational(m, n) * x)
                    coef2[j, i, k, l] += (-1) ** (indicator * (abs(i - j) - 1)) * \
                        sp.exp(sp.Rational(-m, n) * x)

                    # Apply T again
                    if k_human > l_human:
                        k_human, l_human = (T(k_human - 1) % n + 1, T(l_human))
                        k = k_human - 1
                        l = l_human - 1
                        indicator = triple.C((i_human, j_human), (k_human, l_human))
                    else:
                        l_human, k_human = (T(k_human - 1) % n + 1, T(l_human))
                        k = k_human - 1
                        l = l_human - 1
                        indicator = triple.C((i_human, j_human), (k_human, l_human))

                    # counter += 1
                    # print(f"i: {i_human}, j: {j_human}, k: {k_human}, l: {l_human}")
                    # print(f"indicator: {indicator}")
                    # print(f"Here counter is {counter}\n======================")

    return mat2.MatrixTensor2(n, coef1, True) + mat2.MatrixTensor2(n, coef2, True) + \
        int(standard_part) * (1 / (sp.exp(x) - 1)) * mat2.casimir(n) + int(standard_part) * r0