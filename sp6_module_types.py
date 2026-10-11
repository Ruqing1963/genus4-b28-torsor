"""
sp6_module_types.py -- the GL_6(F_2)-conjugacy type (elementary divisors of the F_2[z]-module) of every element of
Sp_6(F_2), with multiplicities, and which types occur in S_8 = O_6^+(2) and in O_6^-(2).

The group is built as in octad_combinatorics.py (closure of S_8 and one tetrad transvection on the octad space).
For every element g the ranks of p(g)^j (p irreducible of degree <= 6, j = 1..6) are computed with numpy from the
permutation table of g on the 64 vectors: the image of a linear map is a subspace, so its rank is log2 of the number
of distinct values in the row.  The invariant has the same format as genus_check.module_invariant, so the output
can be matched against the Ore-algebra computation of the Frobenius module of the genus-4 family.

Output: sp6_module_types.json  = list of {type, label, count, in_S8, in_O6minus, S8_classes}
        (count / 1451520 is the Chebotarev prediction for the frequency of the type among the parameters lambda if
        the monodromy of the family is Sp_6(F_2)).
Usage: python sp6_module_types.py      (about 5 minutes, ~1 GB of memory)
"""
import json
import time
from collections import Counter

import numpy as np

from genus_check import canon_ranks, irreducibles, pdeg, pstr, s8_modules
from octad_combinatorics import OctadModel
from ore_algebra_k6 import divisors_label

DIM = 6
N = 1 << DIM


def rank_of_table(T):
    S = np.sort(T, axis=1)
    distinct = 1 + (np.diff(S, axis=1) != 0).sum(axis=1)
    return np.log2(distinct).round().astype(np.int8)


def main():
    t0 = time.time()
    M = OctadModel(8)
    adj = []
    for i in range(7):
        p = list(range(8))
        p[i], p[i + 1] = p[i + 1], p[i]
        adj.append(tuple(p))
    gens = [M.perm_matrix(p) for p in adj] + [M.transvection(M.coord[0b1111])]
    SP = M.closure(gens)
    P = M.perm_table(SP)  # (G, 64): images of the 64 vectors
    G = len(P)
    print(f"|Sp_6(F_2)| = {G}  ({time.time() - t0:.0f} s)")
    # powers g^i, i = 0..6
    rows = np.arange(G)[:, None]
    powers = [np.tile(np.arange(N, dtype=np.uint8), (G, 1))]
    for i in range(1, DIM + 1):
        powers.append(P[rows, powers[-1]])
    # quadratic forms
    Q0 = np.array(M.Q0, dtype=np.uint8)
    fix_plus = (Q0[P] == Q0[None, :]).all(axis=1)
    minus = None
    for c in range(N):
        vals = np.array([M.Q0[v] ^ M.IP[c][v] for v in range(N)], dtype=np.uint8)
        if (vals == 0).sum() == 28:
            minus = vals
            break
    fix_minus = (minus[P] == minus[None, :]).all(axis=1)
    print(f"fix Q+: {fix_plus.sum()}, fix one Q-: {fix_minus.sum()}")

    # ranks of p(g)^j
    ranks = {}  # p -> array (jmax, G)
    for p in irreducibles(DIM):
        if p == 0b10:
            continue
        jmax = DIM // pdeg(p)
        # table of p(g)
        T = np.zeros((G, N), dtype=np.uint8)
        for i in range(pdeg(p) + 1):
            if p >> i & 1:
                T ^= powers[i]
        R = np.zeros((jmax, G), dtype=np.int8)
        Tj = T
        for j in range(jmax):
            if j:
                Tj = Tj[rows, T]  # p(g)^{j+1} = p(g)^j o p(g)
            R[j] = rank_of_table(Tj)
        ranks[p] = R
        print(f"  p = {pstr(p):16s} done ({time.time() - t0:.0f} s)")

    # assemble invariants in the format of genus_check.module_invariant
    plist = [p for p in irreducibles(DIM) if p != 0b10]
    R_all = {p: ranks[p].T for p in plist}  # (G, jmax)
    types = Counter()
    types_S8 = set()
    types_Om = set()
    inv_of = [None] * G
    for g in range(G):
        inv = []
        for p in plist:
            r = R_all[p][g]
            if r[0] == DIM:
                continue
            inv.append((p, canon_ranks(int(x) for x in r)))
        inv = tuple(inv)
        inv_of[g] = inv
        types[inv] += 1
        if fix_plus[g]:
            types_S8.add(inv)
        if fix_minus[g]:
            types_Om.add(inv)
    print(f"distinct GL_6(F_2)-types in Sp_6(F_2): {len(types)};  in S_8: {len(types_S8)};  in O_6^-(2): {len(types_Om)};"
          f"  ({time.time() - t0:.0f} s)")
    _, inv_to_classes, _ = s8_modules()
    out = []
    for inv, cnt in sorted(types.items(), key=lambda kv: -kv[1]):
        out.append({"type": [[p, list(r)] for p, r in inv], "label": divisors_label(inv), "count": cnt,
                    "fraction": cnt / G, "in_S8": inv in types_S8, "in_O6minus": inv in types_Om,
                    "S8_classes": inv_to_classes.get(inv, [])})
    with open("sp6_module_types.json", "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1)
    print(f"{'count':>8s} {'fraction':>9s}  S8  O6-  S8 classes            elementary divisors")
    for o in out:
        print(f"{o['count']:8d} {o['fraction']:9.5f}  {'+' if o['in_S8'] else '-'}   {'+' if o['in_O6minus'] else '-'}   "
              f"{'/'.join(o['S8_classes']):20s}  {o['label']}")
    n_s8 = sum(o["count"] for o in out if o["in_S8"])
    print(f"elements whose type occurs in S_8: {n_s8} ({n_s8 / G:.4f}); elements lying in a conjugate of S_8: 981820 (0.6764)")
    print(f"total {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
