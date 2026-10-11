"""
genus4_class_stats.py -- the genus-4 "Table 2": distribution of the GL_6(F_2)-type of the q-Frobenius on
K = ker P for all (l5, l3) in F_q^2 (l9 = 1), q = 2, 4, 8, 16, 32, 64, 128 [, 256 with argument 8],
by the Ore-algebra method of ore_algebra_k6.py (no field extensions).

For every parameter the script records the order of Frob_q on K, the elementary divisors of the F_2[z]-module K,
whether the type occurs in S_8 (octad parameter: a Frobenius-stable octad is possible) and, if sp6_module_types.json
is present, whether the type occurs in Sp_6(F_2) at all / in O_6^-(2), together with the Chebotarev prediction
count(type)/|Sp_6(F_2)| for the frequency of the type.
Usage: python genus4_class_stats.py [kmax]       (kmax = 7 by default, i.e. up to q = 128)
"""
import json
import os
import sys
import time
from collections import Counter, defaultdict

from genus_check import GF, canon_ranks, charpoly_from_invariant, pmod, pmul, pstr, s8_modules
from ore_algebra_k6 import DIM, IRR6, Ore, divisors_label, kernel_polynomial, subst


class TGF(GF):
    """GF(2^n) with exp/log tables and a Frobenius table"""

    def __init__(self, n):
        super().__init__(n)
        self.exp = [0] * (2 * self.order)
        self.log = [0] * self.order
        # find a primitive element (for q = 2 it is 1)
        g = 1
        for g in range(2, self.order):
            x, seen = 1, set()
            for i in range(self.order - 1):
                seen.add(x)
                x = GF.mul(self, x, g)
            if len(seen) == self.order - 1:
                break
        x = 1
        for i in range(self.order - 1):
            self.exp[i] = x
            self.log[x] = i
            x = GF.mul(self, x, g)
        for i in range(self.order - 1, 2 * self.order):
            self.exp[i] = self.exp[i - (self.order - 1)]
        self.frob_tab = [[0] * n for _ in range(self.order)]
        for c in range(self.order):
            y = c
            for e in range(n):
                self.frob_tab[c][e] = y
                y = GF.mul(self, y, y)

    def mul(self, a, b):
        if a == 0 or b == 0:
            return 0
        return self.exp[self.log[a] + self.log[b]]

    def inv(self, a):
        return self.exp[(self.order - 1 - self.log[a]) % (self.order - 1)]

    def frob(self, c, e):
        return self.frob_tab[c][e % self.n]


class FastOre(Ore):
    def frob(self, c, e):
        return self.F.frob(c, e)


def frobenius_order_fast(O, P, k, nmax=200):
    """smallest n with tau^{kn} ≡ 1 (mod P on the right): iterate r <- (tau^k r) mod P"""
    r = [1]
    for n in range(1, nmax + 1):
        # left-multiply by tau^k: coefficients get raised to 2^k and shifted
        r = [0] * k + [O.frob(c, k) for c in r]
        r = O.rdivmod(r, P)[1]
        if r == [1]:
            return n
    raise RuntimeError("order not found")


def frobenius_module_fast(O, P, k, order):
    inv = []
    for p in IRR6:
        if p == 0b10:
            continue
        if pmod((1 << order) | 1, p) != 0:  # p must divide z^order + 1
            continue
        ranks = []
        pj = 1
        for j in range(1, DIM + 1):
            pj = pmul(pj, p)
            d = len(O.rgcd(P, subst(pj, k))) - 1
            ranks.append(DIM - d)
            if len(ranks) > 1 and ranks[-1] == ranks[-2]:
                break
        if ranks[0] < DIM:
            inv.append((p, canon_ranks(ranks)))
    return tuple(inv)


def main():
    kmax = int(sys.argv[1]) if len(sys.argv) > 1 else 7
    _, inv_to_classes, _ = s8_modules()
    sp6 = {}
    if os.path.exists("sp6_module_types.json"):
        with open("sp6_module_types.json", encoding="utf-8") as f:
            for o in json.load(f):
                key = tuple((p, tuple(r)) for p, r in o["type"])
                sp6[key] = o
        print(f"loaded sp6_module_types.json: {len(sp6)} types, Chebotarev predictions available")
    else:
        print("sp6_module_types.json not found: run sp6_module_types.py for the Chebotarev predictions")

    grand = {}
    for k in range(1, kmax + 1):
        q = 1 << k
        t0 = time.time()
        F = TGF(k)
        O = FastOre(F)
        stats = Counter()
        orders = defaultdict(Counter)
        for l5 in range(q):
            for l3 in range(q):
                P = kernel_polynomial(F, l5, l3)
                n = frobenius_order_fast(O, P, k)
                inv = frobenius_module_fast(O, P, k, n)
                stats[inv] += 1
                orders[inv][n] += 1
        grand[q] = stats
        total = q * q
        n_octad = sum(c for inv, c in stats.items() if inv in inv_to_classes)
        n_sp6 = sum(c for inv, c in stats.items() if inv in sp6) if sp6 else None
        print("\n" + "=" * 118)
        print(f"q = {q}: {total} parameters (l5, l3), {len(stats)} Frobenius types, {time.time() - t0:.1f} s")
        print(f"octad parameters (type occurs in S_8): {n_octad} = {n_octad / total:.4f}" +
              (f";  type occurs in Sp_6(F_2): {n_sp6} = {n_sp6 / total:.4f}" if sp6 else ""))
        print("=" * 118)
        print(f"{'count':>6s} {'freq':>7s} {'Cheb.':>7s}  {'S_8 class':18s} {'Sp6':4s} {'char poly':26s} {'elementary divisors':40s} orders")
        for inv, c in sorted(stats.items(), key=lambda kv: -kv[1]):
            cls = "/".join(inv_to_classes.get(inv, [])) or "non-octad"
            cheb = f"{sp6[inv]['fraction']:7.4f}" if inv in sp6 else ("   --  " if sp6 else "")
            flag = ("S8" if inv in inv_to_classes else ("O-" if inv in sp6 and sp6[inv]["in_O6minus"] else ("no!" if sp6 else "")))
            print(f"{c:6d} {c / total:7.4f} {cheb}  {cls:18s} {flag:4s} {pstr(charpoly_from_invariant(inv)):26s} "
                  f"{divisors_label(inv):40s} {dict(sorted(orders[inv].items()))}")

    # summary across q
    print("\n" + "=" * 118)
    print("octad proportion by q (Chebotarev prediction for 'type occurs in S_8' and for 'element in a conjugate of S_8' = 0.6764)")
    print("=" * 118)
    pred_type_s8 = sum(o["fraction"] for o in sp6.values() if o["in_S8"]) if sp6 else None
    for q, stats in grand.items():
        total = q * q
        n_octad = sum(c for inv, c in stats.items() if inv in inv_to_classes)
        print(f"  q = {q:4d}: octad {n_octad:6d} / {total:6d} = {n_octad / total:.4f}")
    if pred_type_s8 is not None:
        print(f"  prediction (fraction of Sp_6(F_2) whose GL_6-type occurs in S_8): {pred_type_s8:.4f}")
    # machine-readable counts and polynomial fits in q (separately for odd and even k, as in Part III)
    all_types = sorted({inv for st in grand.values() for inv in st}, key=lambda inv: -max(st.get(inv, 0) for st in grand.values()))
    dump = []
    for inv in all_types:
        dump.append({"type": [[p, list(r)] for p, r in inv], "label": divisors_label(inv),
                     "S8_classes": inv_to_classes.get(inv, []), "charpoly": pstr(charpoly_from_invariant(inv)),
                     "counts": {str(q): st.get(inv, 0) for q, st in grand.items()}})
    with open("genus4_counts.json", "w", encoding="utf-8") as f:
        json.dump(dump, f, indent=1)
    try:
        import sympy as sp
        qs = sp.symbols("q")
        print("\n" + "=" * 118)
        print("polynomial interpolation of the counts (degree <= 2 through the three largest q of the same parity of k;")
        print("'ok' = the remaining points of that parity agree, otherwise the count is not a quadratic in q on that parity)")
        print("=" * 118)
        for inv in all_types:
            cls = "/".join(inv_to_classes.get(inv, [])) or "non-octad"
            line = f"  {cls:18s} {divisors_label(inv):38s}"
            for parity, name in ((1, "odd k"), (0, "even k")):
                pts = [(q, grand[q].get(inv, 0)) for q in sorted(grand) if (q.bit_length() - 1) % 2 == parity]
                if len(pts) < 3:
                    continue
                poly = sp.interpolate(pts[-3:], qs)
                rest = all(poly.subs(qs, q) == c for q, c in pts[:-3])
                line += f" | {name}: {sp.factor(poly)} {'ok' if rest else 'NOT polynomial' if len(pts) > 3 else '(3 pts only)'}"
            print(line)
    except ImportError:
        print("sympy not available: no interpolation")

    # class-by-class comparison at the largest q
    if sp6:
        q = max(grand)
        stats = grand[q]
        total = q * q
        print(f"\nper-type comparison at q = {q}: observed frequency vs Chebotarev prediction (types with prediction > 0.5%)")
        for inv, o in sorted(sp6.items(), key=lambda kv: -kv[1]["fraction"]):
            if o["fraction"] < 0.005:
                continue
            obs = stats.get(inv, 0) / total
            cls = "/".join(o["S8_classes"]) or "non-octad"
            print(f"  {cls:18s} {o['label']:40s} observed {obs:7.4f}  predicted {o['fraction']:7.4f}")
        missing = [inv for inv in stats if inv not in sp6]
        print(f"types observed at q = {q} that do not occur in Sp_6(F_2) at all: {len(missing)}")


if __name__ == "__main__":
    main()
