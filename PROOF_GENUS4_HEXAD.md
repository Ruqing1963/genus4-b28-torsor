# The genus-4 coset-sum theorem: complete proof

Setting. q = 2^k, λ = (λ9, λ5, λ3, λ1) ∈ Fq⁴ with λ9 ≠ 0. For n ≥ 1 and x ∈ F_{2^n},

  𝔮^{(n)}(x) = T_n(x) + 2 Tr_{2^n/2}(λ9x⁹ + λ5x⁵ + λ3x³ + λ1x) ∈ Z/4,   T_n(x) = Tr_{W₂(F_{2^n})/W₂(F₂)}([x]),

and s_n(λ) = Σ_{x∈F_{2^n}} i^{𝔮^{(n)}(x)}. When the field is F_{q^m} we write 𝔮^{(m)} for 𝔮^{(mk)}, as in Parts III–IV. Throughout, 𝔮 ≡ Tr (mod 2), because the first Witt component of the trace of [x] is Tr(x).

Define

  T'(x) = λ9x⁸ + λ5x⁴ + λ3x² + x + (λ3x)^{1/2} + (λ5x)^{1/4} + (λ9x)^{1/8},
  L'(x) = T'(x)⁸ = λ9⁸x⁶⁴ + λ5⁸x³² + λ3⁸x¹⁶ + x⁸ + λ3⁴x⁴ + λ5²x² + λ9x,
  B(x) = B28(x) = λ9⁴x²⁸ + λ5⁴x¹² + λ9²x¹⁰ + λ3⁴x⁴ + λ5²x² + λ9x + 1.

K = ker T' (a 6-dimensional F₂-space, since L' is a separable 2-linearised polynomial of 2-degree 6), C = T'^{-1}(1) = w₀ + K.

Everything below is unconditional. The numerical checks in `b28_identity.py` (Part 3) only guard against algebra slips.

---

## Lemma 1 (polar identity)

For x, y in any F_{2^n} containing Fq:  𝔮(x + y) = 𝔮(x) + 𝔮(y) + 2 Tr(y T'(x)).

Proof. In W₂ of a ring of characteristic 2, [x] + [y] = (x + y, xy) = [x + y] + V[xy], and Tr_{W₂/Z4}(V[a]) = 2Tr(a); hence T(x + y) = T(x) + T(y) + 2Tr(xy). For a Gold monomial, (x+y)^{2^i+1} + x^{2^i+1} + y^{2^i+1} = x^{2^i}y + xy^{2^i}, and Tr(λ x y^{2^i}) = Tr((λx)^{2^{−i}} y). Summing the three monomials and the term 2Tr(xy) gives 2Tr(y T'(x)). ∎

Consequences. K is the radical of the polar form; 𝔮 restricted to K is additive (polar form vanishes), hence takes values in {0, 2}, hence Tr(k) = 0 for every k ∈ K ∩ F_{2^n} (since 𝔮 ≡ Tr mod 2). In particular all points of C ∩ F_{2^n} have the same absolute trace (differences lie in K); this is the unconditional half of the trace-parity theorem of Part IV.

## Theorem 2 (the degree-28 polynomial)

  L'(x) = x⁸ B(x)² + B(x) + 1   in Fq[x].

Proof (one line). x⁸B² = λ9⁸x⁶⁴ + λ5⁸x³² + λ9⁴x²⁸ + λ3⁸x¹⁶ + λ5⁴x¹² + λ9²x¹⁰ + x⁸; adding B cancels x²⁸, x¹², x¹⁰ and leaves λ9⁸x⁶⁴ + λ5⁸x³² + λ3⁸x¹⁶ + x⁸ + λ3⁴x⁴ + λ5²x² + λ9x + 1 = L' + 1. ∎

Corollary 2.1. L' + 1 = B · (x⁸B + 1). The two factors are coprime (a common root would give 1 = 0) and L' + 1 is separable (its derivative is λ9 ≠ 0). Hence C = {B = 0} ⊔ {x⁸B + 1 = 0}, with exactly 28 and 36 points.

Corollary 2.2 (characterisation). Y = x⁸B is a polynomial solution of the Artin–Schreier equation Y² + Y = x⁸(L' + 1); the only other solution is Y + 1, which is not divisible by x⁸. So B is determined by L' alone; the recursion b_d = c_d + b_{(d−8)/2}² (c = coefficients of L' + 1) produces it. The same recursion with x⁴ produces Part III's B = λ5²u⁶ + λ3²u² + λ5u + 1 from the genus-2 L'.

## Theorem 3 (coset-sum theorem, genus 4)

Let m ≥ 1 and ε(w) = +1 if B(w) = 0, ε(w) = −1 otherwise. Then

  s_{2m}(λ) = q^m Σ_{w ∈ C ∩ F_{q^m}} ε(w) i^{𝔮^{(m)}(w)}.

Proof. Write V = F_{q^{2m}}, U = F_{q^m}, 𝔮 = 𝔮^{(2m)}, Tr = Tr_{q^{2m}/2}, φ(x) = x^{q^m}, so Tr_{V/U} = φ + 1.

Step 1 (𝔮 on U). For u ∈ U, T_{2mk}(u) = Tr_{W₂(U)/Z4}([u] + [u]) = Tr_{W₂(U)/Z4}(0, u²) = 2Tr_{q^m/2}(u), and 2Tr_{V/2}(f(u)) = 4Tr_{q^m/2}(f(u)) = 0. So 𝔮(u) = 2Tr_{q^m/2}(u).

Step 2 (cosets of U). For c ∈ V, u ∈ U, by Lemma 1 and Step 1,
  𝔮(c + u) = 𝔮(c) + 2Tr_{q^m/2}(u) + 2Tr(u T'(c)) = 𝔮(c) + 2Tr_{q^m/2}(u(1 + (φ+1)T'(c))),
so Σ_{u∈U} i^{𝔮(c+u)} = q^m i^{𝔮(c)} [ (φ+1)T'(c) = 1 ], and 𝔮 is constant on those cosets. Since T' has coefficients in Fq it commutes with φ, so (φ+1)T'(c) = T'(c + c^{q^m}); the condition is w := c + c^{q^m} ∈ C. As φ + 1 maps V onto U with kernel U, the cosets c + U correspond bijectively to w ∈ U, and
  s_{2m} = q^m Σ_{w ∈ C ∩ U} i^{𝔮(c_w)},  c_w any element with c_w + c_w^{q^m} = w.
It remains to show 𝔮(c) = 𝔮^{(m)}(w) + 2[B(w) ≠ 0].

Step 3 (the Witt value of c). Let N = c·c^{q^m} ∈ U. By transitivity of the Witt trace and φ([c]) = [c^{q^m}],
  T_{2mk}(c) = Tr_{W₂(U)/Z4}([c] + [c^{q^m}]) = Tr_{W₂(U)/Z4}(w, N) = T_{mk}(w) + 2Tr_{q^m/2}(N).
The power sums p_j = c^j + c^{jq^m} of the two conjugates satisfy p_j = w p_{j−1} + N p_{j−2} with p₀ = 0, p₁ = w (characteristic 2), whence
  p₁ = w, p₃ = w³ + Nw, p₅ = w⁵ + Nw³ + N²w, p₉ = w⁹ + Nw⁷ + N²w⁵ + N⁴w
(also p₂ = w², p₄ = w⁴, p₆ = w⁶ + N²w², p₇ = w⁷ + Nw⁵ + N³w, p₈ = w⁸). Therefore
  𝔮(c) = T_{mk}(w) + 2Tr_{q^m/2}(N + λ9p₉ + λ5p₅ + λ3p₃ + λ1p₁) = 𝔮^{(m)}(w) + 2Tr_{q^m/2}(N γ_w),
where, using Tr(N²a) = Tr(N a^{1/2}) and Tr(N⁴a) = Tr(N a^{1/4}),
  γ_w = 1 + λ9w⁷ + (λ9w⁵)^{1/2} + (λ9w)^{1/4} + λ5w³ + (λ5w)^{1/2} + λ3w.
Now
  γ_w² = 1 + λ9²w¹⁴ + λ9w⁵ + (λ9w)^{1/2} + λ5²w⁶ + λ5w + λ3²w²,
  γ_w⁴ = 1 + λ9⁴w²⁸ + λ9²w¹⁰ + λ9w + λ5⁴w¹² + λ5²w² + λ3⁴w⁴ = B(w).
From T'(w) = 1, i.e. λ9w⁸ + λ5w⁴ + λ3w² + w = 1 + (λ3w)^{1/2} + (λ5w)^{1/4} + (λ9w)^{1/8},
  wγ_w = w + λ9w⁸ + λ5w⁴ + λ3w² + (λ9w⁷)^{1/2} + (λ9w⁵)^{1/4} + (λ5w³)^{1/2}
       = 1 + (λ3w)^{1/2} + (λ5w)^{1/4} + (λ9w)^{1/8} + (λ9w⁷)^{1/2} + (λ9w⁵)^{1/4} + (λ5w³)^{1/2},
and squaring twice, (wγ_w)⁴ = 1 + λ3²w² + λ5w + (λ9w)^{1/2} + λ9²w¹⁴ + λ9w⁵ + λ5²w⁶ = γ_w². Since squaring is injective, (wγ_w)² = γ_w, so γ_w ∈ {0, w^{−2}} (w ≠ 0 because T'(0) = 0 ≠ 1), and γ_w = 0 iff γ_w⁴ = B(w) = 0.

If B(w) = 0: 𝔮(c) = 𝔮^{(m)}(w). If B(w) ≠ 0: Nγ_w = c c^{q^m}/w² = v·v^{q^m} with v = c/w; since v + v^{q^m} = 1, v^{q^m} = v + 1, so Nγ_w = v² + v with v ∉ U; the Artin–Schreier polynomial Y² + Y + (v² + v) over U has no root in U, hence Tr_{q^m/2}(v² + v) = 1 and 𝔮(c) = 𝔮^{(m)}(w) + 2. ∎

Remark. The proof is Part IV, Theorem 3.3, verbatim, with two more power sums and one more squaring; the genus-g pattern is γ_w^{2^{n−1}} = B(w), (wγ_w)² = γ_w, for the Gold family with exponents 2^i + 1, i ≤ n (genus 2^{n−1}... here n = 3, genus 4), with B of degree 2^{2n−1} − 2^{n−1} and L' = x^{2^n}B² + B + 1. Only the Gold shape is used (bilinearity of the polar form); nothing about octads.

## Theorem 4 (the coset is a torsor of quadratic forms)

Fix m with C ⊂ F_{q^m} (m = order of Frobenius on K suffices; any multiple works). For w ∈ C define ε(w) as above and, for k ∈ K,
  Q_w(k) = [ε(w + k) ≠ ε(w)] ∈ F₂.
Then:
(i) ω(a, b) := Tr_{q^{2m}/2}(d_b T'(d_a)), with d_a, d_b ∈ F_{q^{2m}} any lifts (d_a + d_a^{q^m} = a, d_b + d_b^{q^m} = b), is a well-defined nondegenerate alternating form on K, independent of m;
(ii) every Q_w is a quadratic form on K with polar form ω;
(iii) Arf(Q_w) = 1 (minus type, 28 zeros) iff ε(w) = +1 iff B(w) = 0; Arf(Q_w) = 0 (plus type, 36 zeros) iff x⁸B(w) + 1 = 0;
(iv) w ↦ Q_w is a bijection from C onto the set of quadratic forms on K with polar form ω, equivariant for the translation action of K (Q_{w+k} = Q_w + ω(k, ·)).

Proof. Let c = c_w, and for k ∈ K ∩ F_{q^m} = K let d = d_k be a lift. By Step 3 of Theorem 3 applied to w and to w + k (with lift c + d), and Lemma 1,
  ε(w + k)ε(w) = i^{𝔮(c + d) − 𝔮(c) − 𝔮^{(m)}(w + k) + 𝔮^{(m)}(w)} = i^{𝔮(d) + 2Tr(dT'(c)) − 𝔮^{(m)}(k) − 2Tr_{q^m/2}(kT'(w))},
and T'(w) = 1. All four terms are even (𝔮(d) ≡ Tr_{V/2}(d) = Tr_{q^m/2}(k) ≡ 𝔮^{(m)}(k) mod 2), so
  Q_w(k) = ½[𝔮(d) − 𝔮^{(m)}(k)] + Tr(dT'(c)) + Tr_{q^m/2}(k)  (mod 2).        (★)
The last term vanishes by the consequence of Lemma 1. The first term is a quadratic form in d (𝔮 is a Z/4-quadratic form on V, and 𝔮^{(m)} on U restricted to K is additive, so ½[𝔮(d) − 𝔮^{(m)}(k)] is quadratic in d with polar form Tr(d′T'(d))), the second is linear in d. Changing the lift d by u ∈ U changes 𝔮(d) by 2Tr_{q^m/2}(u) + 2Tr(uT'(d)) = 2Tr_{q^m/2}(u) + 2Tr_{q^m/2}(u·T'(k)) = 2Tr_{q^m/2}(u), and changes 𝔮^{(m)}(k) not at all and Tr(dT'(c)) by Tr(uT'(c)) = Tr_{q^m/2}(u·T'(w)) = Tr_{q^m/2}(u); the changes cancel, so Q_w(k) is well defined and (★) exhibits it as a quadratic form in k with polar form ω(a, b) = Tr(d_b T'(d_a)); the same computation shows ω is well defined. This proves (ii) and the first half of (i).

ω alternating: 2Tr(dT'(d)) = 𝔮(2d) − 2𝔮(d) ≡ 2𝔮(d) (mod 4), so ω(a, a) = 𝔮(d) mod 2 = Tr_{V/2}(d) = Tr_{q^m/2}(a) = 0 for a ∈ K.
ω nondegenerate: let β(x, y) = Tr(yT'(x)) on V, a symmetric bilinear form with radical K (T' is injective on V/K... more precisely rad β = ker T' ∩ V = K, as the trace form is nondegenerate). Let Ṽ = {d ∈ V : d + d^{q^m} ∈ K} = U + (lifts of K), of dimension mk + 6. Then U ⊂ Ṽ^⊥ (β(u, d) = Tr_{q^m/2}(u·T'(Tr_{V/U}d)) = 0), and dim Ṽ^⊥ = dim V − dim Ṽ + dim(Ṽ ∩ K) = 2mk − (mk + 6) + 6 = mk = dim U, so Ṽ^⊥ = U and rad(β|Ṽ) = Ṽ ∩ U = U. Hence β descends to a nondegenerate form on Ṽ/U ≅ K (via d ↦ d + d^{q^m}), which is ω.
Independence of m: by Theorem 3, ε(w) = [B(w) = 0] does not depend on m, hence neither do Q_w and its polar form.

(iii) Σ_{k∈K} ε(w + k) = ε(w) Σ_k (−1)^{Q_w(k)} = ε(w)(−1)^{Arf Q_w}·2³ (Gauss sum of a nondegenerate quadratic form on F₂⁶). The left side is Σ_{w′∈C} ε(w′) = 28 − 36 = −8 by Corollary 2.1. Hence ε(w)(−1)^{Arf Q_w} = −1.

(iv) Q_{w+k}(k′) = [ε(w + k + k′) ≠ ε(w + k)] = Q_w(k + k′) + Q_w(k) = Q_w(k′) + ω(k, k′). So w ↦ Q_w is K-equivariant between two K-torsors (the quadratic forms with polar form ω form a torsor under K ≅ K* via ω), hence bijective. ∎

Corollary 4.1 (sixty-four-point lemma). If C ⊂ F_{q^m}, then for every λ1, s_{2m}(λ) = ±8 q^m i^{p}, where p = 𝔮^{(m)}(w₀) mod 2 is the common parity of 𝔮^{(m)} on C (𝔮^{(m)} is affine mod 4 on C: 𝔮^{(m)}(w + k) = 𝔮^{(m)}(w) + 𝔮^{(m)}(k) + 2Tr(k) = 𝔮^{(m)}(w) + 2ℓ(k) with ℓ linear). Indeed Σ_C ε(w) i^{𝔮^{(m)}(w)} = ε(w₀) i^{𝔮^{(m)}(w₀)} Σ_k (−1)^{Q_{w₀}(k) + ℓ(k)} = ±8 i^{𝔮^{(m)}(w₀)}. Observed: s₁₄ = ∓1024 i (q = 2, λ = (1,1,1), m = 7) and s₁₈ = ±4096 i (λ = (1,0,0), m = 9).

Corollary 4.2 (what is and is not canonical). The 28 roots of B are the minus-type points of the torsor; the 36 roots of x⁸B + 1 the plus-type points. An octad (eight points of C whose pairwise sums and fourfold sums exhaust K∖0) is the same as a plus-type point P together with a choice making K ≅ F₂⁸_even/𝟙; its 28 duads are the k with Q_P(k) = 1, i.e. P + duads = {B = 0}. The octad itself is not determined by λ (Section 4 of the memo: 0 to 7 Frobenius-stable octads per parameter), and Frobenius permutes an octad iff its class on K is conjugate into S₈ = O₆⁺(2) inside GL(K) (observed proportion → 0.790), whereas B and the torsor exist for every λ.

## What remains open

* The signs of s_{2m}(λ) for the 24 Frobenius types (the analogue of Part IV, Theorems 3.9, 3.17, 3.20): Corollary 4.1 settles the case C ⊂ F_{q^m}; for a proper rational affine subspace C ∩ F_{q^m} the sum is a Gauss sum of Q_{w₀} + ℓ restricted to a subspace, hence 0 or ±2^j, and the sign requires the genus-4 trace parity (value of Tr on C) and the analogue of κ₀.
* The counts of the Frobenius types as functions of q (memo 4.1): observed not polynomial on a parity class of k.
* The Legendre-interval degree to which the family x⁹ corresponds, and its constant.
