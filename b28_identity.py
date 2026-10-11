"""
b28_identity.py -- the canonical degree-28 polynomial of the genus-4 Gold family, for all lambda_9.

Part 1  Symbolic, in F_2[l9, l5, l3][x].  Solving x^8 B^2 + B = L' + 1 coefficient by coefficient (b_d = c_d + b_{(d-8)/2}^2)
        gives the unique polynomial solution
            B_28(x) = l9^4 x^28 + l5^4 x^12 + l9^2 x^10 + l3^4 x^4 + l5^2 x^2 + l9 x + 1,
        and the identity  L'(x) = x^8 B_28(x)^2 + B_28(x) + 1  is checked by expanding (one line).  For l9 = 1 the
        B_28-adic digits of L'+1 are (0, 1, x^8), i.e. L'+1 = B_28 (x^8 B_28 + 1).  The genus-2 case is shown alongside.
Part 2  Numerical: eps = +1 <=> B_28(w) = 0 at q = 16, 32 on random (l5, l3) with l9 = 1, and at q = 4, 8 on random
        (l9, l5, l3) with l9 != 1 (then F_+ = B_28 / l9^4).
Part 3  The ingredients of the proof of the coset-sum theorem, checked numerically on every coset point for a few
        parameters: the power sums p_7 = w^7 + N w^5 + N^3 w, p_9 = w^9 + N w^7 + N^2 w^5 + N^4 w of the two conjugates
        c, c^{q^m} (e_1 = w, e_2 = N), gamma_w^4 = B_28(w), (w gamma_w)^2 = gamma_w, gamma_w in {0, w^{-2}}, and
        eps(w) = +1 iff gamma_w = 0.
Usage: python b28_identity.py
"""
import random

from coset_signs import analyse
from genus_check import GF


# polynomials in F_2[l9, l5, l3][x]: dict xdeg -> set of exponent triples (a9, a5, a3)
def cmul(A, B):
    out = set()
    for a in A:
        for b in B:
            out ^= {(a[0] + b[0], a[1] + b[1], a[2] + b[2])}
    return out


def pmul(P, Q):
    out = {}
    for i, A in P.items():
        for j, B in Q.items():
            C = cmul(A, B)
            if C:
                out[i + j] = out.get(i + j, set()) ^ C
                if not out[i + j]:
                    del out[i + j]
    return out


def padd(P, Q):
    out = {i: set(A) for i, A in P.items()}
    for j, B in Q.items():
        out[j] = out.get(j, set()) ^ B
        if not out[j]:
            del out[j]
    return out


def pdivmod(P, B):
    """B monic in x"""
    db = max(B)
    assert B[db] == {(0, 0, 0)}
    R = {i: set(A) for i, A in P.items()}
    Q = {}
    while R and max(R) >= db:
        d = max(R)
        c = R[d]
        Q[d - db] = c
        R = padd(R, pmul({d - db: c}, B))
    return Q, R


def cstr(A):
    def m(t):
        s = []
        for name, e in zip(("l9", "l5", "l3"), t):
            if e:
                s.append(f"{name}^{e}" if e > 1 else name)
        return " ".join(s) if s else "1"
    terms = sorted(A, reverse=True)
    return " + ".join(m(t) for t in terms) if terms else "0"


def pstr(P):
    parts = []
    for i in sorted(P, reverse=True):
        c = cstr(P[i])
        xs = f"x^{i}" if i > 1 else ("x" if i == 1 else "")
        if c == "1":
            parts.append(xs if xs else "1")
        else:
            parts.append(f"({c}){xs}" if len(P[i]) > 1 else f"{c} {xs}".strip())
    return " + ".join(parts) if parts else "0"


def one(a9=0, a5=0, a3=0):
    return {(a9, a5, a3)}


ONE = {0: one()}
X8 = {8: one()}
# L'(x) + 1 for the genus-4 family, general l9
LP1 = {64: one(8), 32: one(0, 8), 16: one(0, 0, 8), 8: one(), 4: one(0, 0, 4), 2: one(0, 2), 1: one(1), 0: one()}


def solve_B(LP1, shift):
    """the polynomial B with x^shift B^2 + B = LP1 (coefficientwise recursion b_d = c_d + b_{(d-shift)/2}^2)"""
    B = {}
    for d in range(max(LP1) + 1):
        c = set(LP1.get(d, set()))
        if d >= shift and (d - shift) % 2 == 0:
            j = (d - shift) // 2
            c ^= {(2 * a, 2 * b, 2 * e) for (a, b, e) in B.get(j, set())}
        if c:
            B[d] = c
    return B


def part1():
    print("=" * 110)
    print("Part 1: the identity L' = x^8 B_28^2 + B_28 + 1 in F_2[l9, l5, l3][x]")
    print("=" * 110)
    B = solve_B(LP1, 8)
    print(f"solution of x^8 B^2 + B = L' + 1 by the recursion:  B_28 = {pstr(B)}   (degree {max(B)})")
    lhs = padd(pmul(X8, pmul(B, B)), B)  # x^8 B^2 + B, to be compared with L' + 1
    print(f"x^8 B_28^2 + B_28 + 1 == L' : {padd(lhs, LP1) == {}}")
    # the terms that cancel
    print(f"   x^8 B_28^2 = {pstr(pmul(X8, pmul(B, B)))}")
    print(f"   (the monomials x^28, x^12, x^10 cancel between x^8 B_28^2 and B_28)")
    # l9 = 1: division and digits
    B1 = {i: {(0, b, e) for (a, b, e) in A} for i, A in B.items()}
    L1 = {i: {(0, b, e) for (a, b, e) in A} for i, A in LP1.items()}
    Q, R = pdivmod(L1, B1)
    print(f"\nl9 = 1:  (L'+1)/B_28 = {pstr(Q)},  remainder {pstr(R)}")
    digits, T = [], dict(L1)
    while T:
        T, r = pdivmod(T, B1)
        digits.append(r)
    print(f"   B_28-adic digits of L'+1: {[pstr(d) for d in digits]}   =>  L' + 1 = B_28 (x^8 B_28 + 1)")
    # genus 2 for comparison: L' + 1 = l5^4 x^16 + l3^4 x^8 + x^4 + l3^2 x^2 + l5 x + 1, shift 4
    LP1_g2 = {16: one(0, 4), 8: one(0, 0, 4), 4: one(), 2: one(0, 0, 2), 1: one(0, 1), 0: one()}
    B2 = solve_B(LP1_g2, 4)
    lhs2 = padd(pmul({4: one()}, pmul(B2, B2)), B2)
    print(f"\ngenus 2 (same recursion with x^4): B = {pstr(B2)};  x^4 B^2 + B + 1 == L': {padd(lhs2, LP1_g2) == {}}"
          f"   (Part III: B = l5^2 u^6 + l3^2 u^2 + l5 u + 1)")
    return B


def part2():
    print("\n" + "=" * 110)
    print("Part 2: eps = +1 <=> B_28(w) = 0  (F_+ = B_28 / l9^4), random parameters")
    print("=" * 110)
    random.seed(7)
    for k, nsamp, general in ((4, 12, False), (5, 8, False), (2, 6, True), (3, 10, True)):
        q = 1 << k
        Fq = GF(k)
        ok = 0
        for _ in range(nsamp):
            l5, l3 = random.randrange(q), random.randrange(q)
            l9 = random.randrange(2, q) if general else 1
            res = analyse(4, k, l5, l3, verbose=False, lt=l9)
            Fp = res["Fplus"]
            inv4 = Fq.inv(Fq.pow(l9, 4))
            pred = [0] * 29
            pred[28] = 1
            pred[12] = Fq.mul(Fq.pow(l5, 4), inv4)
            pred[10] = Fq.mul(Fq.mul(l9, l9), inv4)
            pred[4] = Fq.mul(Fq.pow(l3, 4), inv4)
            pred[2] = Fq.mul(Fq.mul(l5, l5), inv4)
            pred[1] = Fq.mul(l9, inv4)
            pred[0] = inv4
            same = Fp == pred
            ok += same
            checks = all((res["qf"], res["same_omega"], res["alt"], res["nondeg"], res["arf"], res["intrinsic"]))
            print(f"  q={q:2d}: (l9,l5,l3)=({l9:2d},{l5:2d},{l3:2d}) Frob order {res['m']:2d}, |C_+|={res['sizes'][0]}, "
                  f"F_+ == B_28/l9^4: {same}, torsor checks {checks}")
        print(f"  q={q}: {ok}/{nsamp} agree")


def part3():
    print("\n" + "=" * 110)
    print("Part 3: the proof ingredients on every coset point: power sums, gamma_w^4 = B_28(w), (w gamma_w)^2 = gamma_w")
    print("=" * 110)
    random.seed(11)
    cases = [(1, 1, 1, 1), (1, 1, 0, 0), (2, 3, 1, 2), (3, 5, 2, 7), (3, 1, 6, 3)]  # (k, l9, l5, l3)
    for k, l9, l5, l3 in cases:
        res = analyse(4, k, l5, l3, verbose=False, lt=l9)
        E, fam, C, eps, m = res["E"], res["fam"], res["C"], res["eps"], res["m"]
        L9, L5, L3 = fam.lt, fam.l5, fam.l3
        lift_cols = E.linear_map_matrix(lambda c: c ^ E.frob(c, m * k))
        ok_p, ok_g4, ok_wg, ok_set, ok_eps = True, True, True, True, True
        for w in C:
            c = E.solve(lift_cols, w)
            cq = E.frob(c, m * k)
            N = E.mul(c, cq)
            # power sums
            p7 = E.pow(c, 7) ^ E.pow(cq, 7)
            p9 = E.pow(c, 9) ^ E.pow(cq, 9)
            w2 = E.mul(w, w); w3 = E.mul(w2, w); w5 = E.mul(w3, w2); w7 = E.mul(w5, w2); w9 = E.mul(w7, w2)
            N2 = E.mul(N, N); N3 = E.mul(N2, N); N4 = E.mul(N2, N2)
            if p7 != (w7 ^ E.mul(N, w5) ^ E.mul(N3, w)) or p9 != (w9 ^ E.mul(N, w7) ^ E.mul(N2, w5) ^ E.mul(N4, w)):
                ok_p = False
            # gamma_w
            g = (1 ^ E.mul(L9, w7) ^ E.sqrt(E.mul(L9, w5)) ^ E.sqrt(E.sqrt(E.mul(L9, w))) ^ E.mul(L5, w3)
                 ^ E.sqrt(E.mul(L5, w)) ^ E.mul(L3, w))
            B28w = (E.mul(E.pow(L9, 4), E.pow(w, 28)) ^ E.mul(E.pow(L5, 4), E.pow(w, 12)) ^ E.mul(E.mul(L9, L9), E.pow(w, 10))
                    ^ E.mul(E.pow(L3, 4), E.pow(w, 4)) ^ E.mul(E.mul(L5, L5), w2) ^ E.mul(L9, w) ^ 1)
            if E.pow(g, 4) != B28w:
                ok_g4 = False
            wg = E.mul(w, g)
            if E.mul(wg, wg) != g:
                ok_wg = False
            if g not in (0, E.inv(w2)):
                ok_set = False
            if (eps[w] == 1) != (g == 0):
                ok_eps = False
            # and the Witt difference equals 2 Tr(N gamma_w)
        print(f"  q={1 << k}, (l9,l5,l3)=({l9},{l5},{l3}), Frob order {m}: p7, p9 formulas {ok_p}; gamma^4 = B_28(w) {ok_g4}; "
              f"(w gamma)^2 = gamma {ok_wg}; gamma in {{0, w^-2}} {ok_set}; eps(w)=+1 <=> gamma_w = 0 {ok_eps}")


if __name__ == "__main__":
    part1()
    part2()
    part3()
