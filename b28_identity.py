"""
b28_identity.py -- the canonical degree-28 factor of the genus-4 coset polynomial.

coset_signs.py found, by interpolation over F_2, F_4, F_8, that the points of T'^{-1}(1) with eps(w) = +1 (minus-type
quadratic forms) are the roots of
    B_28(x) = x^28 + l5^4 x^12 + x^10 + l3^4 x^4 + l5^2 x^2 + x + 1        (l9 = 1).
Part 1  exact division in F_2[l5, l3][x]:  L'(x) + 1 = B_28(x) F_36(x), F_36 printed, and the B_28-adic expansion of
        L' + 1 (in genus 2: L' + 1 = B (x^4 B + 1), B = u^6 + l3^2 u^2 + u + 1 for l5 = 1).
Part 2  verification of eps = +1 <=> B_28(w) = 0 at q = 16 and q = 32 on random parameters (big fields up to F_{2^300}).
Usage: python b28_identity.py
"""
import random
import sys

from coset_signs import analyse


# polynomials in F_2[l5, l3][x]: dict xdeg -> set of (a, b) monomials l5^a l3^b
def cmul(A, B):
    out = set()
    for (a, b) in A:
        for (c, d) in B:
            out ^= {(a + c, b + d)}
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
    """B monic in x (leading coefficient {(0,0)})"""
    db = max(B)
    assert B[db] == {(0, 0)}
    R = {i: set(A) for i, A in P.items()}
    Q = {}
    while R and max(R) >= db:
        d = max(R)
        c = R[d]
        Q[d - db] = c
        R = padd(R, pmul({d - db: c}, B))
    return Q, R


def cstr(A):
    def m(a, b):
        s = []
        if a:
            s.append(f"l5^{a}" if a > 1 else "l5")
        if b:
            s.append(f"l3^{b}" if b > 1 else "l3")
        return " ".join(s) if s else "1"
    terms = sorted(A, reverse=True)
    return " + ".join(m(a, b) for a, b in terms) if terms else "0"


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


B28 = {28: {(0, 0)}, 12: {(4, 0)}, 10: {(0, 0)}, 4: {(0, 4)}, 2: {(2, 0)}, 1: {(0, 0)}, 0: {(0, 0)}}
LP1 = {64: {(0, 0)}, 32: {(8, 0)}, 16: {(0, 8)}, 8: {(0, 0)}, 4: {(0, 4)}, 2: {(2, 0)}, 1: {(0, 0)}, 0: {(0, 0)}}
B6 = {6: {(0, 0)}, 2: {(0, 2)}, 1: {(0, 0)}, 0: {(0, 0)}}
LP1_g2 = {16: {(0, 0)}, 8: {(0, 4)}, 4: {(0, 0)}, 2: {(0, 2)}, 1: {(0, 0)}, 0: {(0, 0)}}


def part1():
    print("=" * 110)
    print("Part 1: exact division in F_2[l5, l3][x]")
    print("=" * 110)
    Q, R = pdivmod(LP1_g2, B6)
    print(f"genus 2 (l5 = 1): (L'+1) / B = {pstr(Q)},  remainder {pstr(R)}")
    digits, T = [], dict(LP1_g2)
    while T:
        T, r = pdivmod(T, B6)
        digits.append(r)
    print(f"   B-adic digits of L'+1: {[pstr(d) for d in digits]}")
    Q, R = pdivmod(LP1, B28)
    print(f"\ngenus 4 (l9 = 1): B_28 = {pstr(B28)}")
    print(f"   (L'+1) / B_28 = F_36 = {pstr(Q)}")
    print(f"   remainder = {pstr(R)}   (exact division: {not R})")
    digits, T = [], dict(LP1)
    while T:
        T, r = pdivmod(T, B28)
        digits.append(r)
    print(f"   B_28-adic digits of L'+1 (d_0, d_1, d_2):")
    for i, d in enumerate(digits):
        print(f"      d_{i} = {pstr(d)}")
    # the analogous expression: is F_36 = x^8 B_28 + (something)?
    rest = padd(Q, pmul({8: {(0, 0)}}, B28))
    print(f"   F_36 + x^8 B_28 = {pstr(rest)}")


def part2():
    print("\n" + "=" * 110)
    print("Part 2: eps = +1 <=> B_28(w) = 0 at q = 16, 32 (random parameters)")
    print("=" * 110)
    random.seed(7)
    from genus_check import GF
    for k, nsamp in ((4, 12), (5, 8)):
        q = 1 << k
        Fq = GF(k)
        ok = 0
        for _ in range(nsamp):
            l5, l3 = random.randrange(q), random.randrange(q)
            res = analyse(4, k, l5, l3, verbose=False)
            Fp = res["Fplus"]  # index = degree, coefficients in F_q
            pred = [0] * 29
            pred[28] = 1; pred[12] = Fq.pow(l5, 4); pred[10] = 1; pred[4] = Fq.pow(l3, 4); pred[2] = Fq.mul(l5, l5)
            pred[1] = 1; pred[0] = 1
            same = Fp == pred
            ok += same
            print(f"  q={q}: (l5,l3)=({l5:2d},{l3:2d}) Frob order {res['m']:2d}, |C_+|={res['sizes'][0]}, F_+ == B_28: {same}, "
                  f"torsor checks {all((res['qf'], res['same_omega'], res['alt'], res['nondeg'], res['arf'], res['intrinsic']))}")
        print(f"  q={q}: {ok}/{nsamp} agree")


if __name__ == "__main__":
    part1()
    part2()
