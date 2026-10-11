# octad_lab

Independent exploration (not part of the Legendre-interval series Parts I–IV): the 63-point geometry of an octad, O6+(2) ≅ S8 inside Sp6(F2), and the genus-4 Artin–Schreier family y² + y = λ9x⁹ + λ5x⁵ + λ3x³ + λ1x whose polar operator has a six-dimensional kernel.

* `OCTAD_THEORY.md` — the memo: verified combinatorics, corrections to the initial premise (the x⁷ family is not supersingular and carries no quadratic form; the orbits 28 + 35 belong to O6+(2), not Sp6(F2)), the explicit octad for x⁹ + x⁵ + x³ over F2, the octad obstruction (Frobenius of order 9 for x⁹), and a programme.
* `octad_combinatorics.py` → `octad_combinatorics.txt` — closure, dimension, symplectic and quadratic forms, S8 and Sp6(F2) by explicit group closure (orders 40320 and 1451520), orbits, element orders, the 67.6 % of Sp6(F2) lying in a conjugate of S8; the hexad N = 6 for contrast; the general identity Σ C(N, j) = 2^{N−2} − 1.
* `genus_check.py` → `genus_check.txt` — L-polynomials and Newton polygons of y² + y = x⁷ + … over F2; polar-form bilinearity test; module invariants of the 22 classes of S8 on the octad space; Frobenius module on K = ker T' for the x⁹ family and the octad obstruction table; explicit octad and B8(u) = u⁸ + u⁶ + u + 1.

Requirements: Python 3, numpy. Run from this directory:

```
python octad_combinatorics.py > octad_combinatorics.txt
python genus_check.py > genus_check.txt
```
