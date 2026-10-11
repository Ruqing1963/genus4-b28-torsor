"""
octad_sum_test.py -- the signed coset formula for the genus-4 Gauss sums, read off from the data.

Setting (q = 2^k, lambda = (l9, l5, l3, l1) in F_q^4):
    q^(n)(x) = T_n(x) + 2 Tr_{2^n/2}(l9 x^9 + l5 x^5 + l3 x^3 + l1 x)  in Z/4,   T_n = Witt trace of [x] (Part III/IV),
    s_n = sum_{x in F_{2^n}} i^{q^(n)(x)}.
Part A  the polar identity q(x+y) = q(x) + q(y) + 2 Tr(y T'(x)) with
            T'(x) = l9 x^8 + l5 x^4 + l3 x^2 + x + (l3 x)^{1/2} + (l5 x)^{1/4} + (l9 x)^{1/8}
        is tested exhaustively on F_16 (all lambda in F_16^3) and on F_64 (all lambda in F_8^3), against the variant
        with an extra Tr(x) Tr(y).  This settles the polar operator of the genus-4 family (Programme step 1).
Part B  fibre sums.  For c in F_{q^{2m}} put w = c + c^{q^m} in F_{q^m}.  As in Part IV, the sum of i^{q^(2m)} over the
        coset c + F_{q^m} vanishes unless w lies in T'^{-1}(1), and then equals q^m i^{q^(2m)(c)}; the ratio
            eps(w) = i^{q^(2m)(c) - q^(m)(w)}
        is a sign (+-1) depending only on w.  The script computes eps(w) for every rational point w of the coset
        T'^{-1}(1) cap F_{q^m} without assuming any octad, checks s_{2m} = q^m sum_w eps(w) i^{q^(m)(w)} against the
        brute-force s_{2m}, and then looks at the pattern of eps: how many +1, is the +1 set Frobenius-stable, does it
        coincide with the roots of one of the Frobenius-stable octads (8 points), or with something else.
        Cases: q = 2, (l5,l3) = (1,1), m = 7 (coset in F_128, sum over F_{2^14}); (l5,l3) = (0,0), m = 9 (F_{2^18});
        both for l1 = 0, 1; and the small m = 1..6 for (1,1) where the only rational coset point is w = 1.
Usage: python octad_sum_test.py
"""
import itertools
import time
from collections import Counter

import numpy as np

from genus_check import irreducibles, pdeg
from gf2n import BigGF, f2_rank

I4 = np.array([1, 1j, -1, -1j])


def qform_v(F, x, l9, l5, l3, l1, nconj=None):
    """q^(n)(x) on the array x, n = nconj (default F.N): Witt trace over nconj conjugates + 2 Tr_{2^n/2} f(x)"""
    n = F.N if nconj is None else nconj
    # Witt trace over n conjugates
    s0, s1, c = np.zeros_like(x), np.zeros_like(x), x.copy()
    for _ in range(n):
        s1 ^= F.mul_v(s0, c)
        s0 ^= c
        c = F.mul_v(c, c)
    assert np.all(s0 <= 1) and np.all(s1 <= 1)
    fx = (F.mul_v(F.pow_v(x, 9), l9) ^ F.mul_v(F.pow_v(x, 5), l5) ^ F.mul_v(F.pow_v(x, 3), l3) ^ F.mul_v(x, l1))
    t, y = np.zeros_like(x), fx.copy()
    for _ in range(n):
        t ^= y
        y = F.mul_v(y, y)
    assert np.all(t <= 1)
    return (s0 + 2 * s1 + 2 * t) % 4


def Tprime_v(F, x, l9, l5, l3):
    def sq(a):
        return F.mul_v(a, a)

    def sqrt(a):
        return F.pow_v(a, 1 << (F.N - 1))

    x2 = sq(x); x4 = sq(x2); x8 = sq(x4)
    r1 = sqrt(F.mul_v(x, l3)); r2 = sqrt(sqrt(F.mul_v(x, l5))); r3 = sqrt(sqrt(sqrt(F.mul_v(x, l9))))
    return F.mul_v(x8, l9) ^ F.mul_v(x4, l5) ^ F.mul_v(x2, l3) ^ x ^ r1 ^ r2 ^ r3


def part_A():
    print("=" * 110)
    print("Part A: the polar identity of q(x) = T(x) + 2 Tr(l9 x^9 + l5 x^5 + l3 x^3 + l1 x)")
    print("=" * 110)
    for N, kl in ((4, 4), (6, 3)):
        F = BigGF(N)
        emb = F.embed(kl)
        n = 1 << N
        X, Y = np.meshgrid(np.arange(n), np.arange(n), indexing="ij")
        X, Y = X.ravel(), Y.ravel()
        trv = F.tr_v
        fails_plain, fails_tr, total = 0, 0, 0
        for l9, l5, l3 in itertools.product(emb, repeat=3):
            if l9 == 0:
                continue
            qx = qform_v(F, X, l9, l5, l3, 0)
            qy = qform_v(F, Y, l9, l5, l3, 0)
            qxy = qform_v(F, X ^ Y, l9, l5, l3, 0)
            Tx = Tprime_v(F, X, l9, l5, l3)
            polar = (qxy - qx - qy) % 4
            plain = (2 * trv(F.mul_v(Y, Tx))) % 4
            with_tr = (plain + 2 * (trv(X) & trv(Y))) % 4
            fails_plain += int(np.count_nonzero(polar != plain))
            fails_tr += int(np.count_nonzero(polar != with_tr))
            total += len(X)
        print(f"  F_{1 << N}, lambda in F_{1 << kl}^3 (l9 != 0), {total} evaluations: failures of q(x+y)-q(x)-q(y) = 2Tr(y T'(x)): "
              f"{fails_plain};  of the variant + 2Tr(x)Tr(y): {fails_tr}")


def coset_points(F, m, l9, l5, l3):
    """rational points w of T'^{-1}(1) in F_{2^m} (as a subfield of F = F_{2^{2m}}), by brute force over the subfield"""
    sub = np.array([a for a in range(F.order) if F.in_subfield(a, m)], dtype=np.int64)
    Tp = Tprime_v(F, sub, l9, l5, l3)
    return sub[Tp == 1]


def octads_q2_111(F7):
    """the seven Frobenius-stable octads for lambda = (1,1,1) over F_2 inside F_128 (as sets of ints of F7)"""
    out = []
    for f in irreducibles(7):
        if pdeg(f) != 7:
            continue
        root = next(a for a in range(128) if F7.apply_matrix([F7.pow(a, i) for i in range(8)], f) == 0)
        conj = [F7.pow(root, 1 << i) for i in range(7)]
        if f2_rank(conj) == 7:
            out.append(frozenset(conj + [1]))
    return out


def part_B():
    print("\n" + "=" * 110)
    print("Part B: fibre sums and the sign eps(w) on the rational coset points, q = 2")
    print("=" * 110)
    for (l5, l3), m in (((1, 1), 7), ((0, 0), 9)):
        N = 2 * m
        t0 = time.time()
        F = BigGF(N)
        X = np.arange(F.order, dtype=np.int64)
        W = X ^ F.pow_v(X, 1 << m)  # w = c + c^{2^m}
        rational = coset_points(F, m, 1, l5, l3)
        print(f"\n(l5,l3) = ({l5},{l3}), m = {m}: F_{{2^{N}}}, rational coset points |T'^{{-1}}(1) cap F_{{2^{m}}}| = {len(rational)}")
        # the subfield F_{2^m} as a BigGF of its own is not needed: q^(m)(w) is computed inside F with m conjugates
        for l1 in (0, 1):
            Q = qform_v(F, X, 1, l5, l3, l1)
            vals = I4[Q]
            s2m = vals.sum()
            # fibres
            order = np.argsort(W, kind="stable")
            Ws, Vs = W[order], vals[order]
            bounds = np.flatnonzero(np.diff(Ws)) + 1
            starts = np.concatenate([[0], bounds])
            ends = np.concatenate([bounds, [len(Ws)]])
            fibre = {}
            for s, e in zip(starts, ends):
                fibre[int(Ws[s])] = Vs[s:e].sum()
            rat = set(int(w) for w in rational)
            nonzero_outside = sum(1 for w, v in fibre.items() if w not in rat and abs(v) > 1e-9)
            qm_w = qform_v(F, rational, 1, l5, l3, l1, nconj=m)
            eps = {}
            for w, qv in zip(rational, qm_w):
                ratio = fibre[int(w)] / ((1 << m) * I4[qv])
                assert abs(abs(ratio) - 1) < 1e-9 and abs(ratio.imag) < 1e-9, ratio
                eps[int(w)] = int(round(ratio.real))
            rhs = (1 << m) * sum(eps[int(w)] * I4[qv] for w, qv in zip(rational, qm_w))
            plus = sorted(w for w, e in eps.items() if e == 1)
            minus = sorted(w for w, e in eps.items() if e == -1)
            # Frobenius (x -> x^2) stability of the +1 set
            stable = all((F.mul(w, w) in eps) and eps[F.mul(w, w)] == eps[w] for w in eps)
            print(f"  l1 = {l1}: s_{N} = {s2m.real:+.0f}{s2m.imag:+.0f}i;  fibres outside the coset nonzero: {nonzero_outside};"
                  f"  q^m sum eps(w) i^q(w) = {rhs.real:+.0f}{rhs.imag:+.0f}i  (equal: {abs(s2m - rhs) < 1e-6})")
            print(f"         eps = +1 on {len(plus)} points, -1 on {len(minus)} points;  sign pattern Frobenius-stable: {stable}")
            cnt = Counter((eps[int(w)], int(qv)) for w, qv in zip(rational, qm_w))
            print(f"         distribution of (eps(w), q^(m)(w)): {dict(sorted(cnt.items()))}")
            if (l5, l3) == (1, 1):
                # compare with the seven octads (points of F_128 embedded: the subfield of F_{2^14})
                F7 = BigGF(7)
                octs = octads_q2_111(F7)
                emb7 = F.embed(7)
                plus_set = set(plus)
                for j, oc in enumerate(octs):
                    roots = {emb7[a] for a in oc}
                    inter = len(roots & plus_set)
                    print(f"         octad {j + 1}: roots cap (eps=+1) = {inter} of 8, roots cap (eps=-1) = {8 - inter}")
                # is the +1 set the roots of a Frobenius-stable 8-set with the octad property, for ANY choice?
                if len(plus) == 8:
                    e8 = plus[-1]
                    print(f"         the +1 set has 8 points: sum = {int(np.bitwise_xor.reduce(np.array(plus)))}, "
                          f"rank of e_i + e_8: {f2_rank([w ^ e8 for w in plus[:-1]])}")
        print(f"  ({time.time() - t0:.0f} s)")

    # small m for (1,1): the only rational coset point is w = 1
    print("\n(l5,l3) = (1,1), m = 1..6: coset cap F_{2^m} = {1}; s_{2m} = 2^m eps_m(1) i^{q^(m)(1)}")
    for m in range(1, 7):
        F = BigGF(2 * m)
        X = np.arange(F.order, dtype=np.int64)
        rational = coset_points(F, m, 1, 1, 1)
        assert list(rational) == [1], rational
        line = f"  m = {m}: "
        for l1 in (0, 1):
            s2m = I4[qform_v(F, X, 1, 1, 1, l1)].sum()
            q1 = int(qform_v(F, np.array([1]), 1, 1, 1, l1, nconj=m)[0])
            ratio = s2m / ((1 << m) * I4[q1])
            line += f"l1={l1}: s_{2 * m} = {s2m.real:+.0f}{s2m.imag:+.0f}i, q^(m)(1) = {q1}, eps(1) = {ratio.real:+.0f}{ratio.imag:+.0f}i;  "
        print(line)


if __name__ == "__main__":
    part_A()
    part_B()
