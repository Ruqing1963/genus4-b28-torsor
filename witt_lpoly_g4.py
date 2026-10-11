"""
witt_lpoly_g4.py -- L-polynomials of the genus-4 Witt sums.

s_m(lambda) = sum_{x in F_{q^m}} i^{q^{(m)}_lambda(x)} is minus the trace of Frob^m on the Witt-Artin-Schreier sheaf
L_{psi_4}(x, f_lambda(x)) (rank deg f - 1 = 8, Swan conductor 9 at infinity, pure of weight 1).  With p_m = -s_m the
Newton identities give the elementary symmetric functions e_1..e_8 of the eight Frobenius eigenvalues omega_i, i.e.
    L(T) = prod (1 - omega_i T) = 1 - e_1 T + e_2 T^2 - ... + e_8 T^8   in Z[i][T],
and the predicted p_9, p_10 are compared with the computed s_9, s_10 (this checks the rank).  The normalised
eigenvalues omega_i / sqrt(q) are tested for being roots of unity (Frob^N scalar), and their phases are stored for
the figure.  Output: witt_lpoly_g4.json, witt_lpoly_g4.txt.    Usage: python witt_lpoly_g4.py [kmax=2]
"""
import cmath
import json
import sys
from fractions import Fraction

import numpy as np

from genus_check import GF
from gf2n import BigGF
from octad_sum_test import qform_v
from ore_algebra_k6 import Ore, classify, frobenius_module, frobenius_order, kernel_polynomial
from genus_check import s8_modules

I4 = np.array([1, 1j, -1, -1j])
RANK = 8


def gauss_round(z):
    return complex(round(z.real), round(z.imag))


def newton_e(p, deg):
    """elementary symmetric functions e_1..e_deg from power sums p_1..p_deg (complex, exact Gaussian integers)"""
    e = [1]
    for j in range(1, deg + 1):
        s = 0
        for i in range(1, j + 1):
            s += (-1) ** (i - 1) * e[j - i] * p[i - 1]
        e.append(s / j)
    return e[1:]


def predicted_p(e, p, m):
    """p_m from Newton's identity for m > deg"""
    deg = len(e)
    s = 0
    for i in range(1, deg + 1):
        s += (-1) ** (i - 1) * e[i - 1] * p[m - i - 1]
    return s


def main():
    kmax = int(sys.argv[1]) if len(sys.argv) > 1 else 2
    _, inv_to_classes, _ = s8_modules()
    out = []
    lines = []
    for k in range(1, kmax + 1):
        q = 1 << k
        Fq = GF(k)
        O = Ore(Fq)
        M = 10
        fields = {m: BigGF(m * k) for m in range(1, M + 1)}
        embs = {m: fields[m].embed(k) for m in fields}
        lines.append(f"q = {q}: s_m for m = 1..{M}, L-polynomial of degree {RANK}, check of s_9, s_10")
        for l5 in range(q):
            for l3 in range(q):
                P = kernel_polynomial(Fq, l5, l3)
                n = frobenius_order(O, P, k)
                inv = frobenius_module(O, P, k, n)
                cls = classify(inv, inv_to_classes) or f"non-octad({n})"
                for l1 in range(q):
                    s = []
                    for m in range(1, M + 1):
                        F = fields[m]
                        X = np.arange(F.order, dtype=np.int64)
                        e_ = embs[m]
                        Q = qform_v(F, X, e_[1], e_[l5], e_[l3], e_[l1])
                        s.append(gauss_round(I4[Q].sum()))
                    p = [-v for v in s]
                    e = newton_e(p[:RANK], RANK)
                    e_round = [gauss_round(v) for v in e]
                    integral = all(abs(a - b) < 1e-6 for a, b in zip(e, e_round))
                    pred = [predicted_p(e_round, p, m) for m in (9, 10)]
                    rank_ok = all(abs(pred[i] - p[8 + i]) < 1e-6 for i in range(2))
                    coeffs = [1] + [((-1) ** j) * e_round[j - 1] for j in range(1, RANK + 1)]  # L(T) coefficients
                    roots = np.roots(coeffs[::-1])  # roots of sum coeffs[j] T^j  -> 1/omega
                    omegas = 1 / roots
                    absq = np.abs(omegas) / q ** 0.5
                    phases = sorted(((cmath.phase(w) / (2 * cmath.pi)) % 1 for w in omegas))
                    # roots of unity: smallest N <= 240 with (omega/sqrt q)^N = 1 for all i
                    unit = omegas / q ** 0.5
                    N = next((N for N in range(1, 241) if np.allclose(unit ** N, 1, atol=1e-6)), None)
                    out.append({"q": q, "l5": l5, "l3": l3, "l1": l1, "type": cls, "order": n,
                                "s": [[v.real, v.imag] for v in s], "L": [[c.real, c.imag] for c in coeffs],
                                "integral": integral, "rank8_ok": rank_ok, "abs_ok": bool(np.allclose(absq, 1, atol=1e-6)),
                                "phases": phases, "unity_order": N})
                    lines.append(f"  (l5,l3,l1)=({l5},{l3},{l1}) type {cls:12s} Frob order {n:2d}: e in Z[i] {integral}, "
                                 f"s_9,s_10 predicted {rank_ok}, |omega|=sqrt q {bool(np.allclose(absq, 1, atol=1e-6))}, "
                                 f"(omega/sqrt q)^N = 1 for N = {N};  L(T) = {['%+d%+di' % (c.real, c.imag) for c in coeffs]}")
    with open("witt_lpoly_g4.json", "w", encoding="utf-8") as f:
        json.dump(out, f)
    txt = "\n".join(lines)
    with open("witt_lpoly_g4.txt", "w", encoding="utf-8") as f:
        f.write(txt + "\n")
    print(txt)
    print(f"\nall e in Z[i]: {all(o['integral'] for o in out)};  rank 8 confirmed by s_9, s_10 in all cases: {all(o['rank8_ok'] for o in out)};"
          f"  all |omega| = sqrt q: {all(o['abs_ok'] for o in out)}")


if __name__ == "__main__":
    main()
