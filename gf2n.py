"""
gf2n.py -- GF(2^N) for the octad lab: integer polynomial-basis arithmetic for any N (irreducible found by search),
exp/log tables and vectorised numpy operations when N <= 22, Witt-vector trace of length two, subfield embeddings,
and F_2-linear algebra on N-bit integers.
"""
import numpy as np

from genus_check import IRRED, is_irreducible, pdeg, pmod, pmul


def find_irreducible(N):
    if N in IRRED:
        return IRRED[N]
    for a in range(1, N):
        f = (1 << N) | (1 << a) | 1
        if is_irreducible(f):
            return f
    for a in range(1, N):
        for b in range(a + 1, N):
            for c in range(b + 1, N):
                f = (1 << N) | (1 << a) | (1 << b) | (1 << c) | 1
                if is_irreducible(f):
                    return f
    raise RuntimeError(N)


def prime_factors(n):
    out, p = [], 2
    while p * p <= n:
        if n % p == 0:
            out.append(p)
            while n % p == 0:
                n //= p
        p += 1
    if n > 1:
        out.append(n)
    return out


class BigGF:
    def __init__(self, N, tables=None):
        self.N = N
        self.mod = find_irreducible(N)
        self.order = 1 << N
        self.tables = (N <= 22) if tables is None else tables
        if self.tables:
            self._build_tables()

    # scalar arithmetic ------------------------------------------------------------------------------------------
    def _mul_poly(self, a, b):
        return pmod(pmul(a, b), self.mod)

    def mul(self, a, b):
        if self.tables:
            if a == 0 or b == 0:
                return 0
            return int(self.exp[self.log[a] + self.log[b]])
        return self._mul_poly(a, b)

    def pow(self, a, e):
        r = 1
        while e:
            if e & 1:
                r = self.mul(r, a)
            a = self.mul(a, a)
            e >>= 1
        return r

    def inv(self, a):
        if self.tables:
            return int(self.exp[(self.order - 1 - int(self.log[a])) % (self.order - 1)])
        return self.pow(a, self.order - 2)

    def sqrt(self, a):
        return self.pow(a, 1 << (self.N - 1))

    def frob(self, a, e=1):
        for _ in range(e % self.N):
            a = self.mul(a, a)
        return a

    def tr(self, a):
        s, y = 0, a
        for _ in range(self.N):
            s ^= y
            y = self.mul(y, y)
        assert s in (0, 1)
        return s

    def tr_rel(self, a, d):
        """Tr_{2^N / 2^d}(a)"""
        s, y = 0, a
        for _ in range(self.N // d):
            s ^= y
            y = self.frob(y, d)
        return s

    def tr_n(self, a, n):
        """sum of the first n conjugates a, a^2, ..., a^{2^{n-1}} (the trace of the subfield F_{2^n} when a lies in it)"""
        s, y = 0, a
        for _ in range(n):
            s ^= y
            y = self.mul(y, y)
        return s

    def witt(self, a, nconj=None):
        """Tr_{W_2(F_{2^N})/W_2(F_2)}([a]) in Z/4 = s0 + 2 s1 (Witt addition of the N Teichmueller conjugates);
        with nconj = n the Witt trace of the subfield F_{2^n} when a lies in it"""
        s0, s1, c = 0, 0, a
        for _ in range(self.N if nconj is None else nconj):
            s1 ^= self.mul(s0, c)
            s0 ^= c
            c = self.mul(c, c)
        assert s0 in (0, 1) and s1 in (0, 1)
        return (s0 + 2 * s1) % 4

    # tables and vectorised arithmetic ----------------------------------------------------------------------------
    def _build_tables(self):
        n1 = self.order - 1
        primes = prime_factors(n1)
        g = 2
        while True:
            if all(self.pow_poly(g, n1 // p) != 1 for p in primes):
                break
            g += 1
        self.gen = g
        exp = np.zeros(2 * n1, dtype=np.int64)
        log = np.zeros(self.order, dtype=np.int64)
        x = 1
        for i in range(n1):
            exp[i] = x
            log[x] = i
            x = self._mul_poly(x, g)
        exp[n1:] = exp[:n1]
        self.exp, self.log = exp, log

    def pow_poly(self, a, e):
        r = 1
        while e:
            if e & 1:
                r = self._mul_poly(r, a)
            a = self._mul_poly(a, a)
            e >>= 1
        return r

    def mul_v(self, a, b):
        out = np.zeros(np.broadcast(a, b).shape, dtype=np.int64)
        a = np.broadcast_to(a, out.shape)
        b = np.broadcast_to(b, out.shape)
        m = (a != 0) & (b != 0)
        out[m] = self.exp[self.log[a[m]] + self.log[b[m]]]
        return out

    def pow_v(self, a, e):
        out = np.zeros_like(a)
        m = a != 0
        out[m] = self.exp[(self.log[a[m]] * e) % (self.order - 1)]
        return out

    def tr_v(self, a):
        s, y = np.zeros_like(a), a.copy()
        for _ in range(self.N):
            s ^= y
            y = self.mul_v(y, y)
        return s

    def tr_rel_v(self, a, d):
        s, y = np.zeros_like(a), a.copy()
        for _ in range(self.N // d):
            s ^= y
            y = self.pow_v(y, 1 << d)
        return s

    def witt_v(self, a):
        s0, s1, c = np.zeros_like(a), np.zeros_like(a), a.copy()
        for _ in range(self.N):
            s1 ^= self.mul_v(s0, c)
            s0 ^= c
            c = self.mul_v(c, c)
        assert np.all(s0 <= 1) and np.all(s1 <= 1)
        return (s0 + 2 * s1) % 4

    # subfields ------------------------------------------------------------------------------------------------------
    def embed(self, k):
        """F_{2^k} (as in genus_check.GF(k): ints mod IRRED[k]) -> this field; returns the list of images of 0..2^k-1"""
        assert self.N % k == 0
        if k == 1:
            return [0, 1]
        f = IRRED[k]
        # a root of f in this field: an element of order dividing 2^k - 1 that is a root of f
        e = (self.order - 1) // ((1 << k) - 1)
        y = None
        for z in range(2, 2 + 200):
            cand = self.pow(z, e)
            # evaluate f at cand
            val, p = 0, 1
            for i in range(pdeg(f) + 1):
                if f >> i & 1:
                    val ^= p
                p = self.mul(p, cand)
            if val == 0 and (cand != 1 or k == 1):
                y = cand
                break
        assert y is not None
        images = []
        for a in range(1 << k):
            val, p = 0, 1
            for i in range(k):
                if a >> i & 1:
                    val ^= p
                p = self.mul(p, y)
            images.append(val)
        return images

    def in_subfield(self, a, d):
        return self.frob(a, d) == a

    # F_2-linear algebra -------------------------------------------------------------------------------------------
    def linear_map_matrix(self, func):
        """images of the basis 1, x, ..., x^{N-1} under an F_2-linear func"""
        return [func(1 << i) for i in range(self.N)]

    @staticmethod
    def apply_matrix(cols, v):
        r = 0
        i = 0
        while v:
            if v & 1:
                r ^= cols[i]
            v >>= 1
            i += 1
        return r

    def kernel(self, cols):
        """basis of {v : sum_i v_i cols[i] = 0}"""
        N = self.N
        # row reduce the matrix whose columns are cols: work with augmented (image, preimage) pairs
        pivots = {}  # leading bit -> (image, preimage)
        kernel = []
        for i, c in enumerate(cols):
            img, pre = c, 1 << i
            while img:
                hb = img.bit_length() - 1
                if hb in pivots:
                    pi, pp = pivots[hb]
                    img ^= pi
                    pre ^= pp
                else:
                    pivots[hb] = (img, pre)
                    break
            if img == 0:
                kernel.append(pre)
        return kernel

    def solve(self, cols, target):
        """one v with sum v_i cols[i] = target, or None"""
        pivots = {}
        for i, c in enumerate(cols):
            img, pre = c, 1 << i
            while img:
                hb = img.bit_length() - 1
                if hb in pivots:
                    pi, pp = pivots[hb]
                    img ^= pi
                    pre ^= pp
                else:
                    pivots[hb] = (img, pre)
                    break
        img, pre = target, 0
        while img:
            hb = img.bit_length() - 1
            if hb not in pivots:
                return None
            pi, pp = pivots[hb]
            img ^= pi
            pre ^= pp
        return pre


def f2_rank(vectors):
    rows = [v for v in vectors if v]
    rank = 0
    while rows:
        piv = max(rows)
        hb = piv.bit_length() - 1
        rank += 1
        rows = [r ^ piv if r >> hb & 1 else r for r in rows if r != piv]
        rows = [r for r in rows if r]
    return rank
