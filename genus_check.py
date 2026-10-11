r"""
genus_check.py -- which Artin-Schreier family carries the octad geometry?  (independent of the Legendre repositories)

Part A.  The genus-3 family  y^2 + y = x^7 + l5 x^5 + l3 x^3 + l1 x  over F_2 (all eight coefficient vectors):
         point counts over F_2, F_4, F_8, the L-polynomial, its 2-adic Newton polygon, p-rank and supersingularity.
         Also: is 2 Tr(l7 x^7) a Z/4-quadratic form?  (its polar form must be bilinear; it is not.)
Part B.  The genus-4 family  y^2 + y = l9 x^9 + l5 x^5 + l3 x^3 + l1 x  (all exponents of Gold type 2^i + 1).
         By the derivation of Part III the polar operator is
             T'(x) = l9 x^8 + l5 x^4 + l3 x^2 + x + (l3 x)^{1/2} + (l5 x)^{1/4} + (l9 x)^{1/8},
         T'^8 is a 2-linearised polynomial of 2-degree 6 with 2-associate
             ell(z) = l9^8 z^6 + l5^8 z^5 + l3^8 z^4 + z^3 + l3^4 z^2 + l5^2 z + l9,
         so K = ker T' has F_2-dimension exactly 6 (l9 != 0) and the coset T'^{-1}(1) has 64 points.
         For coefficients in F_2 the Frobenius x -> x^2 acts on K as the companion matrix of ell(z); over F_{2^k} as
         its k-th power.  An octad structure (8 points of the coset whose duads and tetrads exhaust K \ 0) that is
         stable under Frobenius forces Frobenius into a conjugate of S_8 = O_6^+(F_2) inside GL_6(F_2).  The script
         computes, for every (l5, l3) in F_2^2 with l9 = 1 and every k, whether the F_2[z]-module (K, Frob^k) is
         isomorphic to the module of some permutation of S_8 on the 6-dimensional octad space (complete invariant:
         ranks of p(g)^j for all irreducible p), and names the class when it is.
         It also constructs the octad explicitly for (l5, l3) = (1, 1): K = ker Tr_{128/2}, octad = the seven
         conjugates of a normal element of trace 1 together with 1, B_8(u) = (u+1) m_alpha(u); and it checks by brute
         force that for (l5, l3) = (0, 0) (the curve y^2 + y = x^9 + l1 x) Frobenius has order 9 on K, which is not an
         element order of S_8.
Usage: python genus_check.py
"""
import functools
import itertools
from collections import Counter
from fractions import Fraction

from octad_combinatorics import OctadModel

# ----------------------------------------------------------------------------------------------------------------
# F_2[z] on integers (bit i = coefficient of z^i), GF(2^n) in polynomial basis
# ----------------------------------------------------------------------------------------------------------------
IRRED = {1: 0b11, 2: 0b111, 3: 0b1011, 4: 0b10011, 5: 0b100101, 6: 0b1000011, 7: 0b10000011, 8: 0b100011101,
         9: 0b1000010001, 10: 0b10000001001, 11: 0b100000000101, 12: 0b1000001010011}


def pdeg(a):
    return a.bit_length() - 1


def pmul(a, b):
    r = 0
    while b:
        if b & 1:
            r ^= a
        a <<= 1
        b >>= 1
    return r


def pmod(a, m):
    dm = pdeg(m)
    while a and pdeg(a) >= dm:
        a ^= m << (pdeg(a) - dm)
    return a


def pdivmod(a, m):
    q, dm = 0, pdeg(m)
    while a and pdeg(a) >= dm:
        s = pdeg(a) - dm
        q |= 1 << s
        a ^= m << s
    return q, a


def pgcd(a, b):
    while b:
        a, b = b, pmod(a, b)
    return a


def ppow_mod(a, e, m):
    r = 1
    a = pmod(a, m)
    while e:
        if e & 1:
            r = pmod(pmul(r, a), m)
        a = pmod(pmul(a, a), m)
        e >>= 1
    return r


def is_irreducible(f):
    d = pdeg(f)
    if d <= 0:
        return False
    x = 0b10
    xp = x
    for i in range(1, d // 2 + 1):
        xp = pmod(pmul(xp, xp), f)  # x^(2^i) mod f
        if pgcd(f, xp ^ x) != 1:
            return False
    return True


def irreducibles(maxdeg):
    return [f for d in range(1, maxdeg + 1) for f in range(1 << d, 1 << (d + 1)) if is_irreducible(f)]


def factor(f):
    """multiset of irreducible factors (brute force, small degrees)"""
    out = []
    for p in irreducibles(pdeg(f)):
        while pdeg(f) >= pdeg(p):
            q, r = pdivmod(f, p)
            if r:
                break
            out.append(p)
            f = q
        if f == 1:
            break
    return out


def pstr(f):
    if f == 0:
        return "0"
    t = [("1" if i == 0 else "z" if i == 1 else f"z^{i}") for i in range(pdeg(f), -1, -1) if f >> i & 1]
    return "+".join(t)


def porder(f):
    """smallest n with f | z^n + 1 (f(0) = 1)"""
    n = 1
    zn = 0b10
    while True:
        if zn == 1:
            return n
        n += 1
        zn = pmod(pmul(zn, 0b10), f)
        if n > 100000:
            raise RuntimeError


class GF:
    def __init__(self, n):
        self.n, self.m, self.order = n, IRRED[n], 1 << n

    def mul(self, a, b):
        return pmod(pmul(a, b), self.m)

    def pow(self, a, e):
        r = 1
        while e:
            if e & 1:
                r = self.mul(r, a)
            a = self.mul(a, a)
            e >>= 1
        return r

    def tr(self, a):
        s, x = 0, a
        for _ in range(self.n):
            s ^= x
            x = self.mul(x, x)
        return s

    def sqrt(self, a):
        return self.pow(a, 1 << (self.n - 1))

    def inv(self, a):
        return self.pow(a, self.order - 2)


# ----------------------------------------------------------------------------------------------------------------
# matrices over F_2 (rows as ints)
# ----------------------------------------------------------------------------------------------------------------
def mat_mul(A, B, d):
    out = []
    for row in A:
        r = 0
        for i in range(d):
            if row >> i & 1:
                r ^= B[i]
        out.append(r)
    return out


def mat_add(A, B):
    return [a ^ b for a, b in zip(A, B)]


def mat_rank(A):
    rows = [r for r in A if r]
    rank = 0
    while rows:
        piv = max(rows)
        hb = piv.bit_length() - 1
        rank += 1
        rows = [r ^ piv if r >> hb & 1 else r for r in rows if r != piv]
        rows = [r for r in rows if r]
    return rank


def mat_pow(A, e, d):
    R = [1 << i for i in range(d)]
    while e:
        if e & 1:
            R = mat_mul(R, A, d)
        A = mat_mul(A, A, d)
        e >>= 1
    return R


def poly_of_matrix(p, A, d):
    """p(A) for p in F_2[z]"""
    R = [0] * d
    P = [1 << i for i in range(d)]
    for i in range(pdeg(p) + 1):
        if p >> i & 1:
            R = mat_add(R, P)
        P = mat_mul(P, A, d)
    return R


IRR6 = irreducibles(6)


def canon_ranks(ranks):
    """rank sequence of p(A)^j up to and including the first stabilised value, trailing repetitions removed"""
    r = list(ranks)
    while len(r) > 1 and r[-1] == r[-2]:
        r.pop()
    return tuple(r)


def module_invariant(A, d=6):
    """complete isomorphism invariant of the F_2[z]-module (F_2^d, A): ranks of p(A)^j, j = 1..d, all irreducible p
    (canonical form: the strictly decreasing part of the sequence, see canon_ranks)"""
    inv = []
    for p in IRR6:
        PA = poly_of_matrix(p, A, d)
        ranks = []
        P = [1 << i for i in range(d)]
        for j in range(1, d + 1):
            P = mat_mul(P, PA, d)
            ranks.append(mat_rank(P))
            if ranks[-1] == ranks[-2] if len(ranks) > 1 else False:
                break
        if ranks[0] < d:
            inv.append((p, canon_ranks(ranks)))
    return tuple(inv)


def charpoly_from_invariant(inv, d=6):
    """product of p^{e_p} with e_p = (d - rank p(A)^j stabilised) / deg p"""
    f = 1
    for p, ranks in inv:
        e = (d - ranks[-1]) // pdeg(p)
        for _ in range(e):
            f = pmul(f, p)
    return f


def companion(f, d):
    """companion matrix of the monic f of degree d (rows = images of basis vectors 1, z, ..., z^{d-1})"""
    rows = [1 << (i + 1) for i in range(d - 1)]
    rows.append(f & ((1 << d) - 1))  # z^d = f - z^d
    return rows


def partitions(n, maxpart=None):
    if maxpart is None:
        maxpart = n
    if n == 0:
        yield ()
        return
    for k in range(min(n, maxpart), 0, -1):
        for rest in partitions(n - k, k):
            yield (k,) + rest


# ----------------------------------------------------------------------------------------------------------------
def part_A():
    print("=" * 100)
    print("Part A: the genus-3 family y^2 + y = x^7 + l5 x^5 + l3 x^3 + l1 x over F_2")
    print("=" * 100)
    fields = {n: GF(n) for n in range(1, 7)}

    def count(l5, l3, l1, n):
        F = fields[n]
        c = 0
        for x in range(F.order):
            x2 = F.mul(x, x)
            x3 = F.mul(x2, x)
            x5 = F.mul(x3, x2)
            x7 = F.mul(x5, x2)
            f = x7 ^ (x5 if l5 else 0) ^ (x3 if l3 else 0) ^ (x if l1 else 0)
            if F.tr(f) == 0:
                c += 1
        return 1 + 2 * c  # one point at infinity, two points above each x with Tr f(x) = 0

    print(f"{'l5 l3 l1':10s} {'N1 N2 N3':12s} {'L(T) = 1 + a1 T + ... + 8 T^6':44s} {'NP slopes':22s} p-rank supersingular  N4..N6 check")
    for l5, l3, l1 in itertools.product((0, 1), repeat=3):
        Ns = [count(l5, l3, l1, n) for n in range(1, 7)]
        S = [2 ** n + 1 - Ns[n - 1] for n in range(1, 7)]  # power sums of the reciprocal roots
        e1 = S[0]
        e2 = Fraction(e1 * S[0] - S[1], 2)
        e3 = Fraction(e2 * S[0] - e1 * S[1] + S[2], 3)
        a = [1, -e1, e2, -e3, 2 * e2, -4 * e1, 8]
        assert all(x.denominator == 1 for x in map(Fraction, a))
        a = [int(x) for x in a]
        # predicted N4, N5, N6 from the L-polynomial (Newton identities continued)
        e = [1, e1, e2, e3, 2 * e2, 4 * e1, 8]  # elementary symmetric functions of the 6 reciprocal roots
        Spred = list(S[:3])
        for n in range(4, 7):
            s = 0
            for i in range(1, n):
                s += (-1) ** (i - 1) * e[i] * Spred[n - 1 - i]
            s += (-1) ** (n - 1) * n * e[n] if n <= 6 else 0
            Spred.append(int(s))
        ok = all(2 ** n + 1 - Spred[n - 1] == Ns[n - 1] for n in range(4, 7))
        # Newton polygon at 2
        pts = []
        for i, c in enumerate(a):
            if c:
                v = 0
                while c % 2 == 0:
                    c //= 2
                    v += 1
                pts.append((i, v))
        hull = [pts[0]]
        for p in pts[1:]:
            while len(hull) >= 2:
                (x1, y1), (x2, y2) = hull[-2], hull[-1]
                if (y2 - y1) * (p[0] - x1) >= (p[1] - y1) * (x2 - x1):
                    hull.pop()
                else:
                    break
            hull.append(p)
        slopes = []
        for (x1, y1), (x2, y2) in zip(hull, hull[1:]):
            slopes += [Fraction(y2 - y1, x2 - x1)] * (x2 - x1)
        prank = slopes.count(0)
        ss = all(s == Fraction(1, 2) for s in slopes)
        Lstr = " + ".join(f"{c}T^{i}" if i else str(c) for i, c in enumerate(a) if c).replace("+ -", "- ")
        print(f"{l5} {l3} {l1}      {str(Ns[:3]):12s} {Lstr:44s} {str([str(s) for s in slopes]):22s} {prank}      {ss!s:6s}       {ok}")
    print("A supersingular curve needs all slopes 1/2; y^2 + y = x^7 (and every member of the family with coefficients in F_2)")
    print("has p-rank 0 but is NOT supersingular: the Newton polygon has slopes 1/3, 2/3.")

    # polar form of 2 Tr(x^7) is not bilinear
    print("\nIs 2 Tr(l7 x^7) a Z/4-quadratic form?  The polar form P(x,y) = Tr(f(x+y) + f(x) + f(y)) must be additive in x.")
    F = GF(4)
    for name, e in (("x^7", 7), ("x^9", 9), ("x^5", 5), ("x^3", 3)):
        def f(x, e=e):
            return F.pow(x, e)

        def P(x, y):
            return F.tr(f(x ^ y) ^ f(x) ^ f(y))

        fails = sum(1 for x1 in range(16) for x2 in range(16) for y in range(16) if P(x1 ^ x2, y) != P(x1, y) ^ P(x2, y))
        print(f"  f = {name} over F_16: additivity failures of the polar form: {fails} of 4096  -> "
              f"{'quadratic form' if fails == 0 else 'NOT a quadratic form (2-weight of the exponent is 3)'}")


# ----------------------------------------------------------------------------------------------------------------
def s8_modules():
    """for each class of S_8 (partition of 8): the module invariant of its action on the octad space"""
    M = OctadModel(8)
    inv_to_classes = {}
    rows = []
    for lam in partitions(8):
        perm, start = list(range(8)), 0
        for c in lam:
            for i in range(c):
                perm[start + i] = start + (i + 1) % c
            start += c
        imgs = M.perm_matrix(tuple(perm))
        A = list(imgs)  # rows = images of the basis vectors in coordinates
        inv = module_invariant(A)
        name = " ".join(str(c) for c in lam)
        inv_to_classes.setdefault(inv, []).append(name)
        rows.append((name, pstr(charpoly_from_invariant(inv))))
    return M, inv_to_classes, rows


def part_B():
    print("\n" + "=" * 100)
    print("Part B: the genus-4 family y^2 + y = x^9 + l5 x^5 + l3 x^3 + l1 x and the octad obstruction")
    print("=" * 100)
    M, inv_to_classes, rows = s8_modules()
    print("characteristic polynomial of each S_8 class on the 6-dimensional octad space V = F_2^8_even / <1>:")
    for name, cp in rows:
        print(f"   {name:16s} {cp}")
    print(f"({len(rows)} classes of S_8, {len(inv_to_classes)} distinct module types)")

    print("\nFrobenius on K = ker T' for l9 = 1, (l5, l3) in F_2^2: ell(z) = z^6 + l5 z^5 + l3 z^4 + z^3 + l3 z^2 + l5 z + 1")
    for l5, l3 in itertools.product((0, 1), repeat=2):
        ell = (1 << 6) | (l5 << 5) | (l3 << 4) | (1 << 3) | (l3 << 2) | (l5 << 1) | 1
        fac = factor(ell)
        n = porder(ell)
        C = companion(ell, 6)
        print(f"\n  (l5,l3) = ({l5},{l3}):  ell = {pstr(ell)} = {' * '.join('(' + pstr(p) + ')' for p in fac)},"
              f"  K is defined over F_2^{n}, Frobenius x->x^2 has order {n} on K")
        realisable = []
        for k in range(1, n + 1):
            inv = module_invariant(mat_pow(C, k, 6))
            cls = inv_to_classes.get(inv)
            realisable.append((k, cls))
        good = [(k, cls) for k, cls in realisable if cls]
        bad = [k for k, cls in realisable if not cls]
        print(f"    k mod {n} for which Frob^k (= Frobenius of F_{{2^k}}) is the action of an S_8 class: "
              f"{[(k, '/'.join(c)) for k, c in good]}")
        print(f"    k mod {n} with NO Frobenius-stable octad possible (module not realised in S_8): {bad}")

    # explicit octad for (1,1): K = ker Tr_{128/2}
    print("\nExplicit octad for (l5,l3) = (1,1): T'(x)^8 = Tr_{128/2}(x), K = trace-zero hyperplane of F_128.")
    F = GF(7)
    K = [x for x in range(128) if F.tr(x) == 0]
    coset = [x for x in range(128) if F.tr(x) == 1]
    assert len(K) == 64 and len(coset) == 64
    # T'(x) itself, with the square roots: all seven terms
    def Tprime(x):
        return (F.pow(x, 8) ^ F.pow(x, 4) ^ F.mul(x, x) ^ x ^ F.sqrt(x) ^ F.sqrt(F.sqrt(x)) ^ F.sqrt(F.sqrt(F.sqrt(x))))
    assert all(Tprime(x) == F.tr(x) for x in range(128)), "T'(x) = Tr(x) on F_128 for l = (1,1,1)"
    alpha = None
    for a in range(2, 128):
        if F.tr(a) != 1:
            continue
        conj = [F.pow(a, 1 << i) for i in range(7)]
        if mat_rank(conj) == 7:
            alpha = a
            break
    conj = [F.pow(alpha, 1 << i) for i in range(7)]
    octad = conj + [1]
    assert functools.reduce(lambda u, v: u ^ v, octad) == 0
    duads = {a ^ b for a, b in itertools.combinations(octad, 2)}
    tetrads = {a ^ b ^ c ^ d for a, b, c, d in itertools.combinations(octad, 4)}
    triads = {a ^ b ^ c for a, b, c in itertools.combinations(octad, 3)}
    print(f"  alpha = {alpha:#09b} (polynomial basis mod {pstr(IRRED[7]).replace('z', 'x')}), Tr(alpha) = 1, conjugates independent")
    print(f"  octad = seven conjugates of alpha and 1; sum = 0: True;  all eight in the coset T'^{{-1}}(1): {all(F.tr(e) == 1 for e in octad)}")
    print(f"  duads: {len(duads)}, tetrads: {len(tetrads)}, disjoint: {not (duads & tetrads)}, "
          f"duads u tetrads = K \\ 0: {duads | tetrads == set(K) - {0}}")
    print(f"  coset = 8 roots + 56 triads: {set(octad) | triads == set(coset)} ({len(triads)} triads)")
    # B_8(u) = (u+1) m_alpha(u): minimal polynomial of alpha over F_2
    m = 1
    for c in conj:
        m = pmul_gf(m, c, F)
    B8 = pmul_gf(m, 1, F)
    print(f"  m_alpha(u) = {pstr(m)}, B_8(u) = (u+1) m_alpha(u) = {pstr(B8)};  coefficient of u^7 is 0: {not (B8 >> 7 & 1)}")
    print("  Frobenius permutes the octad as a 7-cycle fixing 1: class 7 1 of S_8, as predicted by the module table.")

    # brute-force check of the order-9 obstruction for (0,0): K = ker(x^64 + x^8 + x) in F_512
    print("\nBrute force for (l5,l3) = (0,0), i.e. y^2 + y = x^9 + l1 x:  K = {x in F_512 : x^64 + x^8 + x = 0}")
    F9 = GF(9)
    K9 = [x for x in range(512) if F9.pow(x, 64) ^ F9.pow(x, 8) ^ x == 0]
    orders = Counter()
    for x in K9:
        if x == 0:
            continue
        y, o = F9.mul(x, x), 1
        while y != x:
            y = F9.mul(y, y)
            o += 1
        orders[o] += 1
    print(f"  |K| = {len(K9)};  Frobenius orbit lengths on K \\ 0: {dict(sorted(orders.items()))};  "
          f"order of Frobenius on K = {max(orders)} (not an element order of S_8: S_8 has no 9-cycle)")
    print("  Hence over F_{2^k} with 3 does not divide k the curve y^2 + y = x^9 + l1 x admits no Frobenius-stable octad;")
    print("  the module table above says the same for 3 | k, 9 does not divide k, and allows the trivial class 1^8 for 9 | k.")


def pmul_gf(poly, root, F):
    """multiply the polynomial poly (coefficients in F, stored as list) by (u + root); here all coefficients end in F_2
    because we only use it for the full conjugate set, so we keep an integer encoding with coefficients in F"""
    # poly is a dict degree -> coefficient (in F) encoded as int bits only when coefficients are in F_2; use list form
    if isinstance(poly, int):
        coeffs = [(poly >> i) & 1 for i in range(poly.bit_length())]
    else:
        coeffs = poly
    new = [0] * (len(coeffs) + 1)
    for i, c in enumerate(coeffs):
        new[i + 1] ^= c
        new[i] ^= F.mul(c, root)
    if all(c in (0, 1) for c in new):
        return sum(c << i for i, c in enumerate(new))
    return new


if __name__ == "__main__":
    part_A()
    part_B()
