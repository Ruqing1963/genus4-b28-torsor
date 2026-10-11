# octad_lab

Independent exploration (not part of the Legendre-interval series Parts I–IV): the 63-point geometry of an octad, O6+(2) ≅ S8 inside Sp6(F2), and the genus-4 Artin–Schreier family y² + y = λ9x⁹ + λ5x⁵ + λ3x³ + λ1x whose polar operator has a six-dimensional kernel.

* `OCTAD_THEORY.md` — the memo: verified combinatorics, corrections to the initial premise (the x⁷ family is not supersingular and carries no quadratic form; the orbits 28 + 35 belong to O6+(2), not Sp6(F2)), the explicit octad for x⁹ + x⁵ + x³ over F2, the octad obstruction (Frobenius of order 9 for x⁹), and a programme.
* `octad_combinatorics.py` → `octad_combinatorics.txt` — closure, dimension, symplectic and quadratic forms, S8 and Sp6(F2) by explicit group closure (orders 40320 and 1451520), orbits, element orders, the 67.6 % of Sp6(F2) lying in a conjugate of S8; the hexad N = 6 for contrast; the general identity Σ C(N, j) = 2^{N−2} − 1.
* `genus_check.py` → `genus_check.txt` — L-polynomials and Newton polygons of y² + y = x⁷ + … over F2; polar-form bilinearity test; module invariants of the 22 classes of S8 on the octad space; Frobenius module on K = ker T' for the x⁹ family and the octad obstruction table; explicit octad and B8(u) = u⁸ + u⁶ + u + 1.

* `ore_algebra_k6.py` → `ore_algebra_k6.txt` — right division and right gcd in Fq{τ}; the GL6(F2)-type of Frob_q on K = ker P from τ-degrees of right gcds, no field extensions; self-test against the companion matrices (q = 2, 4, 8) and the sixteen parameters of F4.
* `sp6_module_types.py` → `sp6_module_types.json/.txt` — GL6(F2)-types of all 1451520 elements of Sp6(F2) with multiplicities (24 types; 20 in S8, 22 in O6−(2)); Chebotarev predictions.
* `genus4_class_stats.py` → `genus4_class_stats.txt`, `genus4_counts.json` — the genus-4 Frobenius table for all (λ5, λ3) ∈ Fq², q = 2 … 512: type counts, S8 class names, Frobenius orders, comparison with the Sp6(F2) prediction, polynomial fits.

* `gf2n.py` — GF(2^N) for any N (irreducible by search, tables for N ≤ 22, vectorised numpy ops, Witt trace, subfield embeddings, F2-linear algebra on N-bit ints).
* `octad_sum_test.py` → `octad_sum_test.txt` — the polar identity of the genus-4 Witt form verified exhaustively (F16, F64); fibre sums over F_{2^14}, F_{2^18}: the sign ε(w) on the 64 coset points is +1 on 28 and −1 on 36 points, not 8/56.
* `find_octic_identity.py` → `find_octic_identity.txt` — octic search (negative): the seven Frobenius-stable octads of λ = (1,1,1), their B8 and B8-adic expansions; factorisation of L'+1 over F2 (no degree-8 factor for λ = (0,0)); all Frobenius-stable octads over F4, F8 with their coefficients; rational coset points = 2^{dim ker(Frob−1)}.
* `coset_signs.py` → `coset_signs_g2.txt`, `coset_signs_g4.txt` — ε from one lift per point; the quadratic-form torsor (Q_P quadratic, canonical ω = Tr(d_b T'(d_a)), Arf rule) for all λ at q = 2, 4, 8 in genus 2 (recovers B) and genus 4; interpolation of the canonical factor F₊ = B28.
* `b28_identity.py` → `b28_identity.txt` — exact division L' + 1 = B28 (x⁸ B28 + 1) in F2[λ5, λ3][x]; verification of ε = +1 ⟺ B28(w) = 0 at q = 16, 32.

Requirements: Python 3, numpy, sympy. Run from this directory:

```
python octad_combinatorics.py > octad_combinatorics.txt
python genus_check.py > genus_check.txt
python ore_algebra_k6.py > ore_algebra_k6.txt
python sp6_module_types.py > sp6_module_types.txt      # 2 minutes, 1 GB
python genus4_class_stats.py 9 > genus4_class_stats.txt # q up to 512, 6 minutes
python octad_sum_test.py > octad_sum_test.txt
python find_octic_identity.py 3 > find_octic_identity.txt
python coset_signs.py 2 3 > coset_signs_g2.txt
python coset_signs.py 4 3 > coset_signs_g4.txt             # 5 minutes
python b28_identity.py > b28_identity.txt                  # 10 minutes (q = 32 in F_{2^300})
```

Main result so far: the genus-4 coset T'^{-1}(1) is canonically a torsor of quadratic forms on K = F2⁶; its 28 minus-type points are the roots of B28(x) = x²⁸ + λ5⁴x¹² + x¹⁰ + λ3⁴x⁴ + λ5²x² + x + 1, with L' = x⁸B28² + B28 + 1, the exact analogue of the hexad polynomial of Part III.
