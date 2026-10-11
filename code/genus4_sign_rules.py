r"""
genus4_sign_rules.py -- the Gauss sums s_{2m}(lambda) on the 24 Frobenius types, from the coset-sum theorem.

By the coset-sum theorem, s_{2m} = q^m sum_{w in C_m} eps(w) i^{q^{(m)}(w)}, C_m = C cap F_{q^m} = w_r + K_m
(K_m = K cap F_{q^m}, dim f_m; C_m is never empty), eps(w) = +1 iff B_28(w) = 0.  By the torsor theorem and the
polar identity, q^{(m)}(w_r + k) = q^{(m)}(w_r) + 2 l(k) with l linear on K_m (l depends on lambda_1 through
Tr_{q^m/2}(lambda_1 k)), and eps(w_r + k) = eps(w_r) (-1)^{Q(k)}, Q = Q_{w_r}; hence

      s_{2m} / q^m = eps(w_r) i^{q^{(m)}(w_r)} sum_{k in K_m} (-1)^{Q(k) + l(k)} ,

the Gauss sum of the quadratic form Q + l restricted to the subspace K_m.  With r_m = dim rad(omega | K_m):
      s_{2m} = 0  if Q + l is not identically 0 on rad(omega|K_m),  otherwise |s_{2m}| = 2^{(f_m + r_m)/2} q^m,
and the set of lambda_1 in F_q with s_{2m} != 0 is an affine subspace (of F_q as an F_2-space) or empty or all.
For m = d t with d = gcd(m, n): C_m = C_d and q^{(m)} = t q^{(d)} on C_d (multiplicity lemma), so only d | n matters.

The script computes, for all (l5, l3) in F_q^2 (l9 = 1), all lambda_1 in F_q, all d | n: f_d, r_d, the values
s_{2d}/q^d (exact, from the coset data; no sums over F_{q^{2d}}), checks the structure above, and tabulates by
Frobenius type: the observed value sets, the number of lambda_1 with s != 0, and whether the sign (the +- in
+-2^{(f+r)/2} q^d i^p) is constant over the parameters of the type.  For q = 2 the values are compared with the
brute-force sums of octad_sum_test.py (s_14, s_18) and with s_2, ..., s_12.
Usage: python genus4_sign_rules.py [kmax=3]
"""
import sys
from collections import Counter, defaultdict

from coset_signs import analyse
from genus_check import GF, s8_modules
from gf2n import f2_rank
from ore_algebra_k6 import Ore, classify, frobenius_module, frobenius_order, kernel_polynomial

I4 = [1, 1j, -1, -1j]


def divisors(n):
    return [d for d in range(1, n + 1) if n % d == 0]


def fmt(z):
    z = complex(round(z.real, 6), round(z.imag, 6))
    if z == 0:
        return "0"
    if z.imag == 0:
        return f"{int(z.real):+d}"
    if z.real == 0:
        return f"{int(z.imag):+d}i"
    return f"{int(z.real):+d}{int(z.imag):+d}i"


def main():
    kmax = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    _, inv_to_classes, _ = s8_modules()
    for k in range(1, kmax + 1):
        q = 1 << k
        Fq = GF(k)
        O = Ore(Fq)
        print("=" * 118)
        print(f"q = {q}: s_{{2d}}/q^d on the Frobenius types, d | n, all lambda_1 in F_q")
        print("=" * 118)
        table = defaultdict(lambda: defaultdict(Counter))  # type -> d -> Counter of (f, r, nonzero count, value set)
        agg = {}
        struct_ok = True
        for l5 in range(q):
            for l3 in range(q):
                P = kernel_polynomial(Fq, l5, l3)
                n = frobenius_order(O, P, k)
                inv = frobenius_module(O, P, k, n)
                cls = classify(inv, inv_to_classes) or "non-octad(" + ("9" if n % 9 == 0 else str(n)) + ")"
                res = analyse(4, k, l5, l3, verbose=False)
                E, fam, C, K, eps, omega = res["E"], res["fam"], res["C"], res["K"], res["eps"], res["omega"]
                emb = fam.emb
                for d in divisors(n):
                    Cd = [w for w in C if E.frob(w, d * k) == w]
                    Kd = [a for a in K if E.frob(a, d * k) == a]
                    fd = len(Kd).bit_length() - 1
                    # radical of omega restricted to K_d
                    rad = [a for a in Kd if all(omega[(a, b)] == 0 for b in Kd)]
                    rd = len(rad).bit_length() - 1
                    # q^{(d)}_0(w) for w in C_d (lambda_1 = 0), and traces for the lambda_1 term
                    q0 = {w: fam.q_sub(w, d * k) for w in Cd}
                    trw = {w: E.tr_n(w, d * k) for w in Cd}  # Tr_{q^d/2}(w)
                    vals = {}
                    for l1 in range(q):
                        L1 = emb[l1]
                        s = 0
                        for w in Cd:
                            qv = (q0[w] + 2 * E.tr_n(E.mul(L1, w), d * k)) % 4
                            s += eps[w] * I4[qv]
                        vals[l1] = complex(round(s.real), round(s.imag))
                    nz = [l1 for l1 in range(q) if vals[l1] != 0]
                    mags = {abs(vals[l1]) for l1 in nz}
                    expected = 2 ** ((fd + rd) / 2)
                    ok = all(abs(mg - expected) < 1e-6 for mg in mags) and (len(nz) == 0 or (len(nz) & (len(nz) - 1)) == 0)
                    struct_ok &= ok
                    valset = tuple(sorted({fmt(vals[l1]) for l1 in range(q)}))
                    # the lambda_1-sum, analogue of W' of Parts II-IV: W_d = sum_{l1} s_{2d} / q^{d+1}
                    W = sum(vals.values()) / q
                    table[cls][d][(fd, rd, len(nz), valset, fmt(W))] += 1
                    agg.setdefault(cls, {}).setdefault(str(d), {"f": fd, "r": rd, "nz": Counter(), "mag": expected,
                                                                "W": Counter(), "n": 0})
                    a = agg[cls][str(d)]
                    a["nz"][len(nz)] += 1
                    a["W"][fmt(W)] += 1
                    a["n"] += 1
                    if k == 1 and (l5, l3) in ((1, 1), (0, 0)) and d == n:
                        print(f"  check vs brute force: (l5,l3)=({l5},{l3}), m={n}: s_{{2m}}/2^m = "
                              f"{[fmt(vals[l1]) for l1 in range(q)]} (octad_sum_test: -8i,+8i / +8i,-8i)")
        print(f"  structure |s_{{2d}}| in {{0, 2^((f+r)/2) q^d}} and #{{l1: s != 0}} a power of 2: {struct_ok}")
        for cls in sorted(table, key=lambda c: (c.startswith("non"), c)):
            print(f"\n  type {cls}:")
            for d in sorted(table[cls]):
                for (fd, rd, nnz, valset, W), cnt in sorted(table[cls][d].items()):
                    print(f"     d={d:2d}  f={fd} r={rd}  #l1 with s!=0: {nnz}/{q}  values of s_{{2d}}/q^d: {{{', '.join(valset)}}}"
                          f"  W_d = sum_l1 s_{{2d}}/q^(d+1) = {W}  x{cnt} parameters")
        import json
        import os
        DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
        with open(os.path.join(DATA, f"genus4_sign_rules_q{q}.json"), "w", encoding="utf-8") as f:
            json.dump({cls: {d: {"f": a["f"], "r": a["r"], "mag": a["mag"], "n": a["n"], "nz": dict(a["nz"]), "W": dict(a["W"])}
                             for d, a in dd.items()} for cls, dd in agg.items()}, f, indent=1)


if __name__ == "__main__":
    main()
