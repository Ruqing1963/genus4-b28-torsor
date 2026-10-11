"""
find_octic_identity.py -- is there an octic B_8(u; lambda) behind the genus-4 kernel polynomial
    L'(x) = x^64 + l5^8 x^32 + l3^8 x^16 + x^8 + l3^4 x^4 + l5^2 x^2 + x      (l9 = 1),
as B_lambda(u) = l5^2 u^6 + l3^2 u^2 + l5 u + 1 is behind L' = x^4 B^2 + B + 1 in genus 2 (Part III)?

Part 0  genus-2 sanity check of the method: the B-adic expansion of (L'+1)/B is x^4 B + 1.
Part 1  lambda = (1,1,1) over F_2: all Frobenius-stable octads (normal elements of F_128 of trace 1, up to conjugacy),
        their octics B_8 = (u+1) m_alpha(u), the B_8-adic expansion of the triad polynomial T_56 = (L'+1)/B_8, and the
        test of the pattern L' = x^32 B_8^4 + x^b B_8^2 + B_8 + 1.
Part 2  no canonical octic: the factorisation of L'+1 over F_2 for the four lambda in F_2^2 (sympy) -- for
        lambda = (0,0) the degrees are 1, 9, ..., 9, so no degree-8 factor, hence no B_8(u; lambda) in F_2[lambda][u].
Part 3  the Frobenius-stable octads of every octad parameter over F_4 and F_8: construction of K and of the coset
        T'^{-1}(1) by F_2-linear algebra in F_{q^n}, n = order of Frobenius on K; enumeration of the octads from the
        Frobenius orbits on the coset (cycle type = S_8 class); the elementary symmetric functions of the eight roots.
Part 4  rational points of the coset: #{w in F_q : L'(w) = 1} = 2^{dim ker(Frob_q - 1)|K} for every lambda, q <= 64
        (the coset is an affine K-set on which Frobenius always has a fixed point).
Usage: python find_octic_identity.py [kmax for part 3, default 3]
"""
import itertools
import sys
from collections import Counter
from math import comb

from genus_check import GF, IRRED, factor, irreducibles, pdeg, pdivmod, pmul, pstr, s8_modules
from gf2n import BigGF, f2_rank
from ore_algebra_k6 import Ore, classify, frobenius_module, frobenius_order, kernel_polynomial, subst


def expand_base(T, B):
    """B-adic digits d_0, d_1, ... with T = sum d_i B^i, deg d_i < deg B (polynomials over F_2 as ints)"""
    digits = []
    while T:
        T, r = pdivmod(T, B)
        digits.append(r)
    return digits


def poly_mul_list(coeffs, F):
    """(u + r_1)...(u + r_n) over a BigGF field F; returns list of coefficients, index = degree"""
    P = [1]
    for r in coeffs:
        Q = [0] * (len(P) + 1)
        for i, c in enumerate(P):
            Q[i + 1] ^= c
            Q[i] ^= F.mul(c, r)
        P = Q
    return P


def part0():
    print("=" * 110)
    print("Part 0: genus 2, lambda = (1,1) over F_2: B = u^6+u^2+u+1, L' = x^16+x^8+x^4+x^2+x")
    print("=" * 110)
    B = 0b1000111
    Lp = (1 << 16) | (1 << 8) | (1 << 4) | (1 << 2) | (1 << 1)
    q, r = pdivmod(Lp ^ 1, B)
    print(f"  (L'+1)/B = {pstr(q)}, remainder {r};  B-adic digits of (L'+1)/B: "
          f"{[pstr(d) for d in expand_base(q, B)]}   (expected [1, x^4]: (L'+1)/B = x^4 B + 1)")


def part1():
    print("\n" + "=" * 110)
    print("Part 1: lambda = (1,1,1) over F_2: the seven Frobenius-stable octads and their octics")
    print("=" * 110)
    F = BigGF(7)
    Lp = (1 << 64) | (1 << 32) | (1 << 16) | (1 << 8) | (1 << 4) | (1 << 2) | (1 << 1)
    septics = [f for f in irreducibles(7) if pdeg(f) == 7]
    found = []
    for m in septics:
        # a root of m in F_128: brute force
        root = next(a for a in range(128) if F.apply_matrix([F.pow(a, i) for i in range(8)], m) == 0)
        conj = [F.pow(root, 1 << i) for i in range(7)]
        normal = f2_rank(conj) == 7
        if not normal:
            continue
        B8 = pmul(m, 0b11)
        T56, rem = pdivmod(Lp ^ 1, B8)
        assert rem == 0
        digits = expand_base(T56, B8)
        found.append((m, B8, digits))
        # the user's pattern: L' = x^32 B^4 + x^b B^2 + B + 1  <=>  T56 = x^32 B^3 + x^b B + 1
        patt = (len(digits) == 4 and digits[0] == 1 and digits[2] == 0 and digits[3] == (1 << 32)
                and bin(digits[1]).count("1") == 1)
        print(f"  m_alpha = {pstr(m):22s} B_8 = {pstr(B8):28s} T_56 = (L'+1)/B_8 has B_8-adic digits (d_0, d_1, ...):")
        for i, d in enumerate(digits):
            print(f"      d_{i} = {pstr(d)}")
        print(f"      pattern x^32 B^3 + x^b B + 1: {patt}")
    print(f"  {len(septics)} irreducible septics, {len(found)} with normal roots = {len(found)} Frobenius-stable octads "
          f"(all roots of a normal septic have trace 1 automatically)")
    print("  Conclusion: none of the seven octics satisfies an identity of the genus-2 shape; the digits are not monomials.")


def part2():
    print("\n" + "=" * 110)
    print("Part 2: factorisation of the coset polynomial L'(x) + 1 over F_2 for lambda in F_2^2 (sympy)")
    print("=" * 110)
    import sympy as sp
    x = sp.symbols("x")
    for l5, l3 in itertools.product((0, 1), repeat=2):
        Lp1 = x ** 64 + l5 * x ** 32 + l3 * x ** 16 + x ** 8 + l3 * x ** 4 + l5 * x ** 2 + x + 1
        _, facs = sp.factor_list(Lp1, modulus=2)
        degs = sorted((sp.degree(f, x) for f, e in facs for _ in range(e)), reverse=True)
        print(f"  (l5,l3)=({l5},{l3}): factor degrees {degs}   (= Frobenius orbit sizes on the 64-point coset)")
    print("  For (0,0) the degrees are 1 and 9: no product of factors has degree 8, so no octic B_8(u; l5, l3) with")
    print("  polynomial coefficients can divide L'+1 for all lambda.  (For (1,1) the degree-1 factor is x+1 = u+1 and the")
    print("  seven septics of Part 1 are among the degree-7 factors.)")


def part3(kmax):
    print("\n" + "=" * 110)
    print("Part 3: Frobenius-stable octads for every octad parameter over F_q, q = 4, 8")
    print("=" * 110)
    _, inv_to_classes, _ = s8_modules()
    for k in range(2, kmax + 1):
        q = 1 << k
        Fq = GF(k)
        O = Ore(Fq)
        print(f"\n--- q = {q} (F_q = F_2[w]/({pstr(IRRED[k]).replace('z', 'w')}), elements written as integers 0..{q - 1}) ---")
        summary = Counter()
        for l5 in range(q):
            for l3 in range(q):
                P = kernel_polynomial(Fq, l5, l3)
                n = frobenius_order(O, P, k)
                inv = frobenius_module(O, P, k, n)
                cls = classify(inv, inv_to_classes)
                if cls is None:
                    continue
                N = k * n
                F = BigGF(N)
                emb = F.embed(k)
                back = {v: i for i, v in enumerate(emb)}
                L5, L3 = emb[l5], emb[l3]

                def Lprime(x):
                    x2 = F.mul(x, x); x4 = F.mul(x2, x2); x8 = F.mul(x4, x4); x16 = F.mul(x8, x8)
                    x32 = F.mul(x16, x16); x64 = F.mul(x32, x32)
                    return (x64 ^ F.mul(F.pow(L5, 8), x32) ^ F.mul(F.pow(L3, 8), x16) ^ x8 ^ F.mul(F.pow(L3, 4), x4)
                            ^ F.mul(F.mul(L5, L5), x2) ^ x)

                cols = F.linear_map_matrix(Lprime)
                K = F.kernel(cols)
                assert len(K) == 6, (q, l5, l3, len(K))
                w0 = F.solve(cols, 1)
                assert w0 is not None
                coset = []
                for mask in range(64):
                    v = w0
                    for i in range(6):
                        if mask >> i & 1:
                            v ^= K[i]
                    coset.append(v)
                assert all(Lprime(w) == 1 for w in coset)
                # Frobenius orbits on the coset
                seen, orbits = set(), []
                for w in coset:
                    if w in seen:
                        continue
                    orb, y = [], w
                    while y not in orb:
                        orb.append(y)
                        y = F.frob(y, k)
                    seen |= set(orb)
                    orbits.append(orb)
                by_len = {}
                for orb in orbits:
                    by_len.setdefault(len(orb), []).append(orb)
                orbit_sizes = sorted((len(o) for o in orbits), reverse=True)
                # enumerate octads for each cycle type of the class
                octads = []
                for cname in cls.split("/"):
                    ctype = Counter(int(c) for c in cname.split())
                    if any(len(by_len.get(L, [])) < mult for L, mult in ctype.items()):
                        continue
                    ncomb = 1
                    for L, mult in ctype.items():
                        ncomb *= comb(len(by_len[L]), mult)
                    if ncomb > 200000:
                        print(f"  (l5,l3)=({l5},{l3}) class {cname}: {ncomb} combinations, skipped")
                        continue
                    choices = [itertools.combinations(by_len[L], mult) for L in ctype]
                    for pick in itertools.product(*choices):
                        pts = [w for group in pick for orb in group for w in orb]
                        if len(pts) != 8:
                            continue
                        s = 0
                        for w in pts:
                            s ^= w
                        if s != 0:
                            continue
                        e8 = pts[-1]
                        if f2_rank([w ^ e8 for w in pts[:-1]]) != 6:
                            continue
                        octads.append((cname, pts))
                summary[(cls, len(octads))] += 1
                coeff_strs = []
                for cname, pts in octads:
                    B8 = poly_mul_list(pts, F)
                    assert all(c in back for c in B8), "coefficients not in F_q?"
                    b = [back[c] for c in B8]  # b[i] = coefficient of u^i
                    assert b[8] == 1 and b[7] == 0
                    coeff_strs.append("(" + ",".join(str(b[i]) for i in range(6, -1, -1)) + ")")
                print(f"  (l5,l3)=({l5},{l3}) class {cls:10s} Frob order {n:2d}, coset orbit sizes {orbit_sizes}, "
                      f"Frobenius-stable octads: {len(octads)}")
                if coeff_strs:
                    print(f"        B_8 coefficients (b_6,b_5,b_4,b_3,b_2,b_1,b_0): {' '.join(coeff_strs)}")
        print(f"  summary (class, number of octads) -> number of parameters: {dict(sorted(summary.items(), key=str))}")


def part4():
    print("\n" + "=" * 110)
    print("Part 4: rational points of the coset versus the Frobenius module, all lambda, q <= 64")
    print("=" * 110)
    for k in range(1, 7):
        q = 1 << k
        Fq = GF(k)
        O = Ore(Fq)
        ok, mn, dist = True, 99, Counter()
        for l5 in range(q):
            for l3 in range(q):
                l58, l38, l34, l52 = Fq.pow(l5, 8), Fq.pow(l3, 8), Fq.pow(l3, 4), Fq.mul(l5, l5)
                cnt = 0
                for w in range(q):
                    w2 = Fq.mul(w, w); w4 = Fq.mul(w2, w2); w8 = Fq.mul(w4, w4); w16 = Fq.mul(w8, w8)
                    w32 = Fq.mul(w16, w16); w64 = Fq.mul(w32, w32)
                    if (w64 ^ Fq.mul(l58, w32) ^ Fq.mul(l38, w16) ^ w8 ^ Fq.mul(l34, w4) ^ Fq.mul(l52, w2) ^ w) == 1:
                        cnt += 1
                P = kernel_polynomial(Fq, l5, l3)
                d = len(O.rgcd(P, subst(0b11, k))) - 1  # dim ker(Frob_q - 1) on K
                ok &= cnt == (1 << d)
                mn = min(mn, cnt)
                dist[cnt] += 1
        print(f"  q = {q:3d}: #(coset cap F_q) = 2^dim ker(Frob-1) for all lambda: {ok};  minimum {mn};  "
              f"distribution {dict(sorted(dist.items()))}")


if __name__ == "__main__":
    kmax = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    part0()
    part1()
    part2()
    part3(kmax)
    part4()
