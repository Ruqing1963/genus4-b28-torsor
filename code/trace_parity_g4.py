r"""
trace_parity_g4.py -- the trace parity theorem in genus 4 (unconditional), checked on every rational point.

Theorem (trace parity).  For every m >= 1 and every w in C cap F_{q^m}, C = T'^{-1}(1):
      Tr_{q^m/2}(w)  =  m k + f_m   (mod 2),       f_m = dim_{F_2} (K cap F_{q^m}).

Proof.  Fix m and put U = F_{q^m}, N = mk, K_m = K cap U.  For mu in U let q'(x) = q^{(m)}_lambda(x) + 2 Tr_{U/F_2}(mu x);
this is a Z/4-quadratic form on U with the same polar form 2 Tr(y T'(x)) (the linear term drops out), so its radical
is K_m, and q' is additive on K_m with values in {0, 2}.  Choose mu with Tr_{U/F_2}(mu k) = q'_0(k)/2 on K_m (possible:
k -> Tr(mu k) runs through all linear forms of U, hence of K_m); then q' vanishes on K_m and, by Brown's theorem on
Gauss sums of Z/4-forms, s' = sum_{x in U} i^{q'(x)} has |s'|^2 = 2^{N + f_m} and s' != 0.
Shift by w in C cap U: q'(x + w) = q'(x) + q'(w) + 2 Tr(w T'(x)) = q'(x) + q'(w) + 2 Tr(x T'(w)) = q'(x) + q'(w) + 2 Tr(x)
(symmetry of the polar form, T'(w) = 1), and i^{2 Tr(x)} = (-1)^{q'(x)} = i^{-2 q'(x)}, so
      s' = i^{q'(w)} conj(s').
Hence s'^2 = i^{q'(w)} 2^{N + f_m}.  Writing s' = a + b i with a, b integers: if q'(w) is even, s' is real or purely
imaginary and N + f_m is even; if q'(w) is odd, a = +-b and N + f_m is odd.  So q'(w) = N + f_m (mod 2), and
q'(w) = Tr_{U/F_2}(w) (mod 2) because the first Witt component of the Witt trace is the trace.  qed

(The same argument makes Theorem 3.5 of Part IV unconditional: the auxiliary linear term mu may be taken in F_{q^m}
rather than in F_q.)   Consequently all points of C cap F_{q^m} have trace mk + f_m, and K cap F_{q^m} has trace 0.

The script verifies the congruence for every lambda = (1, l5, l3) over F_q, q = 2, 4, 8, every m = 1..n (n = order of
Frobenius on K; C cap F_{q^m} = C cap F_{q^gcd(m,n)}), every rational point, and reports the result by Frobenius type,
with f_m and the trace value.  Usage: python trace_parity_g4.py [kmax=3]
"""
import sys
from collections import Counter, defaultdict

from coset_signs import analyse
from genus_check import s8_modules
from ore_algebra_k6 import classify, frobenius_module, frobenius_order, kernel_polynomial, Ore
from genus_check import GF


def main():
    kmax = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    _, inv_to_classes, _ = s8_modules()
    for k in range(1, kmax + 1):
        q = 1 << k
        Fq = GF(k)
        O = Ore(Fq)
        print("=" * 112)
        print(f"q = {q}: Tr_{{q^m/2}}(w) = mk + f_m (mod 2) for all w in C cap F_{{q^m}}, all m, all (l5, l3)")
        print("=" * 112)
        total_pts, bad = 0, 0
        per_type = defaultdict(Counter)
        for l5 in range(q):
            for l3 in range(q):
                P = kernel_polynomial(Fq, l5, l3)
                n = frobenius_order(O, P, k)
                inv = frobenius_module(O, P, k, n)
                cls = classify(inv, inv_to_classes) or "non-octad:" + str(hash(inv) % 1000)
                res = analyse(4, k, l5, l3, verbose=False)
                E, C, K = res["E"], res["C"], res["K"]
                for m in range(1, n + 1):
                    Cm = [w for w in C if E.frob(w, m * k) == w]
                    Km = [a for a in K if E.frob(a, m * k) == a]
                    fm = len(Km).bit_length() - 1
                    assert 1 << fm == len(Km) and len(Cm) == len(Km)
                    traces = {E.tr_n(w, m * k) for w in Cm}
                    total_pts += len(Cm)
                    ok = traces == {(m * k + fm) % 2}
                    if not ok:
                        bad += len(Cm)
                    per_type[cls][(m, fm, tuple(sorted(traces)), ok)] += 1
        print(f"  rational points checked: {total_pts};  violations: {bad}")
        for cls in sorted(per_type, key=lambda c: (c.startswith("non"), c)):
            items = per_type[cls]
            # compress: for each m, (f_m, trace)
            summary = {}
            for (m, fm, tr, ok), cnt in sorted(items.items()):
                summary.setdefault(m, set()).add((fm, tr[0] if len(tr) == 1 else tr, ok))
            line = "; ".join(f"m={m}: f={sorted({s[0] for s in v})} Tr={sorted({s[1] for s in v})}{'' if all(s[2] for s in v) else ' VIOLATION'}"
                             for m, v in summary.items())
            print(f"  type {cls:12s} {line}")


if __name__ == "__main__":
    main()
