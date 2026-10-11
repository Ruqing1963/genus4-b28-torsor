r"""
coset_signs.py -- the canonical sign eps(w) on the coset T'^{-1}(1), the quadratic-form torsor, and the canonical
factor F_+(x; lambda) of the coset polynomial, for the Gold families of genus 2 (x^5, x^3) and genus 4 (x^9, x^5, x^3).

Facts used (polar identity verified in octad_sum_test.py; coset argument of Part IV):
  for w in C = T'^{-1}(1) cap F_{q^m} and any c in F_{q^{2m}} with c + c^{q^m} = w, the number
      eps(w) = i^{ q^{(2m)}(c) - q^{(m)}(w) }
  is +-1 and depends only on w (not on c, not on lambda_1); the fibre sum over c + F_{q^m} equals q^m eps(w) i^{q^(m)(w)}.
So eps(w) is computed here from ONE lift c per point, by F_2-linear algebra in E = F_{q^{2m}}, m = order of Frobenius on K.
For each lambda the script
  * builds K = ker T' and C in F_{q^m} (inside E), computes eps on all 2^{2n} points of C (2n = dim K),
  * checks: eps is Frobenius-stable; |C_+| = 2^{2n-1} - 2^{n-1} (minus-type count), |C_-| = 2^{2n-1} + 2^{n-1};
  * defines Q_P(k) = [eps(P + k) != eps(P)] on K for each P in C and checks that Q_P is a quadratic form whose polar
    form omega(a, b) = Q_P(a+b) + Q_P(a) + Q_P(b) is alternating, nondegenerate and independent of P, with
    Arf(Q_P) = 1 iff eps(P) = +1; and that omega(a, b) = Tr_E(d_b T'(d_a)) for lifts d_a + d_a^{q^m} = a;
  * forms F_+(x) = prod_{w in C_+} (x - w) in F_q[x] (degree 6 in genus 2, where it must be B_lambda; degree 28 in genus 4)
    and F_- = (L'+1)/F_+;
and finally interpolates every coefficient of F_+ as a reduced polynomial in (l5, l3) over F_q (degrees < q), to find
the closed form of the canonical factor (genus 2: B = l5^2 u^6 + l3^2 u^2 + l5 u + 1 must come out).
Usage: python coset_signs.py [genus=4] [kmax=3]
"""
import itertools
import sys
import time
from collections import Counter

from genus_check import GF, IRRED, pstr
from gf2n import BigGF, f2_rank
from ore_algebra_k6 import Ore, frobenius_order

I4 = [1, 1j, -1, -1j]


def setup(genus, Fq, l5, l3, lt=1):
    """Ore kernel polynomial P (tau-coefficients) and the number n with dim K = 2n; lt = l9 (genus 4)"""
    if genus == 4:
        P = [lt, Fq.mul(l5, l5), Fq.pow(l3, 4), 1, Fq.pow(l3, 8), Fq.pow(l5, 8), Fq.pow(lt, 8)]
        return P, 3
    # genus 2: T'(x) = l5 x^4 + l3 x^2 + x + (l3 x)^{1/2} + (l5 x)^{1/4};  T'^4 = l5^4 x^16 + l3^4 x^8 + x^4 + l3^2 x^2 + l5 x
    P = [l5, Fq.mul(l3, l3), 1, Fq.pow(l3, 4), Fq.pow(l5, 4)]
    return P, 2


class Family:
    def __init__(self, genus, E, k, l5, l3, lt=1):
        self.genus, self.E, self.k = genus, E, k
        emb = E.embed(k)
        self.emb = emb
        self.back = {v: i for i, v in enumerate(emb)}
        self.l5, self.l3, self.lt = emb[l5], emb[l3], emb[lt]

    def sq(self, a):
        return self.E.mul(a, a)

    def Tp(self, x):
        E, l5, l3 = self.E, self.l5, self.l3
        x2 = self.sq(x); x4 = self.sq(x2)
        r1 = E.sqrt(E.mul(l3, x)); r2 = E.sqrt(E.sqrt(E.mul(l5, x)))
        if self.genus == 2:
            return E.mul(l5, x4) ^ E.mul(l3, x2) ^ x ^ r1 ^ r2
        x8 = self.sq(x4)
        r3 = E.sqrt(E.sqrt(E.sqrt(E.mul(self.lt, x))))
        return E.mul(self.lt, x8) ^ E.mul(l5, x4) ^ E.mul(l3, x2) ^ x ^ r1 ^ r2 ^ r3

    def f(self, x):
        E = self.E
        x2 = self.sq(x); x3 = E.mul(x2, x); x5 = E.mul(x3, x2)
        v = E.mul(self.l5, x5) ^ E.mul(self.l3, x3)
        if self.genus == 4:
            x4 = self.sq(x2)
            v ^= E.mul(self.lt, E.mul(x5, x4))
        return v

    def q_big(self, c):
        """q^{(N)}(c) on E"""
        return (self.E.witt(c) + 2 * self.E.tr(self.f(c))) % 4

    def q_sub(self, w, nconj):
        """q^{(nconj)}(w) for w in the subfield of degree nconj"""
        return (self.E.witt(w, nconj) + 2 * self.E.tr_n(self.f(w), nconj)) % 4


def analyse(genus, k, l5, l3, verbose=True, lt=1):
    Fq = GF(k)
    O = Ore(Fq)
    P, n = setup(genus, Fq, l5, l3, lt)
    m = frobenius_order(O, P, k)
    N = 2 * m * k
    E = BigGF(N, tables=False)
    fam = Family(genus, E, k, l5, l3, lt)
    Lp = lambda x: E.frob(fam.Tp(x), n)  # L' = T'^{2^n}
    cols = E.linear_map_matrix(Lp)
    K = E.kernel(cols)
    assert len(K) == 2 * n, (len(K), 2 * n)
    assert all(E.in_subfield(v, m * k) for v in K)
    w0 = E.solve(cols, 1)
    assert w0 is not None
    # move w0 into F_{q^m}: the coset has a point there (find by trying w0 + K)
    C = None
    for mask in range(1 << (2 * n)):
        v = w0
        for i in range(2 * n):
            if mask >> i & 1:
                v ^= K[i]
        if E.in_subfield(v, m * k):
            w0 = v
            break
    assert E.in_subfield(w0, m * k), "coset has no point in F_{q^m}?"
    C = []
    for mask in range(1 << (2 * n)):
        v = w0
        for i in range(2 * n):
            if mask >> i & 1:
                v ^= K[i]
        C.append(v)
    assert all(fam.Tp(w) == 1 for w in C)
    # lifts c + c^{q^m} = w
    lift_cols = E.linear_map_matrix(lambda c: c ^ E.frob(c, m * k))
    eps = {}
    for w in C:
        c = E.solve(lift_cols, w)
        assert c is not None
        e = (fam.q_big(c) - fam.q_sub(w, m * k)) % 4
        assert e in (0, 2), e
        eps[w] = 1 if e == 0 else -1
    Cp = [w for w in C if eps[w] == 1]
    Cm = [w for w in C if eps[w] == -1]
    stable = all(eps[E.frob(w, k)] == eps[w] for w in C)
    # quadratic forms Q_P
    Kset = set()
    for mask in range(1 << (2 * n)):
        v = 0
        for i in range(2 * n):
            if mask >> i & 1:
                v ^= K[i]
        Kset.add(v)
    Klist = sorted(Kset)
    omega_ref = None
    qf_ok, same_omega, arf_ok = True, True, True
    for Pt in C:
        Q = {a: int(eps[Pt ^ a] != eps[Pt]) for a in Klist}
        om = {(a, b): Q[a ^ b] ^ Q[a] ^ Q[b] for a in Klist for b in Klist}
        # bilinear?
        for a in Klist[:8]:
            for b in Klist:
                for c in Klist[:8]:
                    if om[(a ^ c, b)] != om[(a, b)] ^ om[(c, b)]:
                        qf_ok = False
        if omega_ref is None:
            omega_ref = om
        elif om != omega_ref:
            same_omega = False
        zeros = sum(1 for a in Klist if Q[a] == 0)
        arf = 0 if zeros == (1 << (2 * n - 1)) + (1 << (n - 1)) else 1
        if (arf == 1) != (eps[Pt] == 1):
            arf_ok = False
    alt = all(omega_ref[(a, a)] == 0 for a in Klist)
    nondeg = all(any(omega_ref[(a, b)] for b in Klist) for a in Klist if a)
    # intrinsic omega: Tr_E(d_b T'(d_a))
    lifts = {a: E.solve(lift_cols, a) for a in Klist}
    Tl = {a: fam.Tp(lifts[a]) for a in Klist}
    sub_a = Klist if N < 60 else Klist[:12]  # the full check is slow in the big fields of q = 8
    intrinsic_ok = all(omega_ref[(a, b)] == E.tr(E.mul(lifts[b], Tl[a])) for a in sub_a for b in Klist)
    # the canonical factor
    def poly_from_roots(roots):
        Pc = [1]
        for r in roots:
            Qc = [0] * (len(Pc) + 1)
            for i, cc in enumerate(Pc):
                Qc[i + 1] ^= cc
                Qc[i] ^= E.mul(cc, r)
            Pc = Qc
        return Pc
    Fp = poly_from_roots(Cp)
    assert all(cc in fam.back for cc in Fp)
    Fp_q = [fam.back[cc] for cc in Fp]
    res = dict(m=m, sizes=(len(Cp), len(Cm)), stable=stable, qf=qf_ok, same_omega=same_omega, alt=alt, nondeg=nondeg,
               arf=arf_ok, intrinsic=intrinsic_ok, Fplus=Fp_q, E=E, fam=fam, C=C, eps=eps, K=Klist)
    if verbose:
        print(f"  (l5,l3)=({l5},{l3}): Frob order {m:2d}, |C_+|={len(Cp)}, |C_-|={len(Cm)}, eps Frobenius-stable {stable}, "
              f"Q_P quadratic {qf_ok}, polar form independent of P {same_omega}, alternating {alt}, nondegenerate {nondeg}, "
              f"Arf(Q_P)=1 <=> eps(P)=+1 {arf_ok}, omega = Tr(d_b T'(d_a)) {intrinsic_ok}")
        print(f"        F_+ coefficients (x^0 .. x^{len(Fp_q) - 1}): {Fp_q}")
    return res


def solve_fq(Fq, M, rhs):
    """solve M c = rhs over F_q (M list of rows), Gaussian elimination"""
    nrow, ncol = len(M), len(M[0])
    A = [row[:] + [r] for row, r in zip(M, rhs)]
    piv_cols = []
    r = 0
    for cidx in range(ncol):
        p = next((i for i in range(r, nrow) if A[i][cidx]), None)
        if p is None:
            continue
        A[r], A[p] = A[p], A[r]
        inv = Fq.inv(A[r][cidx])
        A[r] = [Fq.mul(inv, v) for v in A[r]]
        for i in range(nrow):
            if i != r and A[i][cidx]:
                fct = A[i][cidx]
                A[i] = [v ^ Fq.mul(fct, u) for v, u in zip(A[i], A[r])]
        piv_cols.append(cidx)
        r += 1
    sol = [0] * ncol
    for i, cidx in enumerate(piv_cols):
        sol[cidx] = A[i][ncol]
    return sol


def interpolate(genus, k, results):
    """reduced polynomial in (l5, l3) over F_q of every coefficient of F_+"""
    Fq = GF(k)
    q = 1 << k
    pts = sorted(results)
    deg = len(results[pts[0]]["Fplus"]) - 1
    monos = [(a, b) for a in range(q) for b in range(q)]
    M = [[Fq.mul(Fq.pow(u, a), Fq.pow(v, b)) for (a, b) in monos] for (u, v) in pts]
    print(f"\n  reduced polynomials of the coefficients of F_+ over F_{q} (monomials l5^a l3^b with 0 <= a, b < q,"
          f" coefficients in F_q as ints;\n  an exponent a >= (q-1)/2 may also be read as a - (q-1), e.g. l5^{q - 3} = l5^(-2) on l5 != 0):")

    def mono_str(a, b):
        s = []
        if a:
            s.append(f"l5^{a}")
        if b:
            s.append(f"l3^{b}")
        return " ".join(s) if s else "1"

    for j in range(deg, -1, -1):
        rhs = [results[p]["Fplus"][j] for p in pts]
        sol = solve_fq(Fq, M, rhs)
        terms = [f"{'' if c == 1 else str(c) + '*'}{mono_str(a, b)}" for (a, b), c in zip(monos, sol) if c]
        if terms:
            print(f"    x^{j:2d}: {' + '.join(terms)}")


def main():
    genus = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    kmax = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    print("=" * 112)
    print(f"genus {genus}: canonical signs on the coset, quadratic-form torsor, canonical factor F_+")
    print("=" * 112)
    for k in range(1, kmax + 1):
        q = 1 << k
        print(f"\n--- q = {q} (F_q = F_2[w]/({pstr(IRRED[k]).replace('z', 'w')})) ---")
        t0 = time.time()
        results = {}
        agg = Counter()
        for l5 in range(q):
            for l3 in range(q):
                if genus == 2 and l5 == 0:
                    continue  # genus 2 needs l5 != 0 (dim K = 4)
                res = analyse(genus, k, l5, l3)
                results[(l5, l3)] = res
                agg[(res["sizes"], res["stable"], res["qf"], res["same_omega"], res["alt"], res["nondeg"], res["arf"],
                     res["intrinsic"])] += 1
        print(f"  summary ((|C+|,|C-|), stable, quadratic, same omega, alternating, nondegenerate, Arf rule, intrinsic omega) -> count:")
        for key, cnt in agg.items():
            print(f"     {key}: {cnt}")
        interpolate(genus, k, results)
        print(f"  ({time.time() - t0:.0f} s)")


if __name__ == "__main__":
    main()
