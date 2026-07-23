# The Correlation-Damping Theorem for the Krasikov Hierarchy

**Formal write-up of link 3 of the limit-theorem proof program (Obs 321–325).
Status per component marked explicitly. 2026-07.**

> **CORRECTION NOTICE (Obs 328).** Theorem A below is VACUOUS in the
> infinite-unrolling limit: every infinite equation-tree path eventually
> accumulates arbitrarily many feed edges, so M_g(m) → v(m) as the unrolling
> depth grows (verified numerically: T_g/T_{g−1} → 1). The finite-depth
> version retains content only when supplemented with leaf-oscillation
> bounds. The correct functional for the damping mechanism is the
> **second-moment (L²) influence propagation** of the linearized Perron
> operator, not first-moment mass; Lemmas 1–2 are unaffected and remain the
> combinatorial core. Two further consequences: (i) the exact flow identity
> φ̄ = 1 − λ⁻²/ρ (see below) is a genuine one-line theorem (backbone is a
> permutation, so its share of total flow is exactly λ⁻²/ρ) and fully
> explains the measured feed fraction 0.7049; (ii) the three-decimal match
> between that number and the measured covariance-decay ratio 0.70 (Obs 324)
> is therefore possibly coincidental until the L² mechanism is worked out —
> downgraded from "confirmed mechanism" to "suggestive". This notice
> supersedes the affected claims; the file is kept in original form below
> for the record.

Throughout: classes are odd residues m ≡ 2 (mod 3) taken mod 3^k; the Krasikov
system I_k has, per class, the fixed-point (Perron) equation

    v(m) = ρ⁻¹ [ λ⁻² v(4m) + φ-term(m) ],

where the feed term is λ^(α−2)·v̄(r₁(m)) for D1 (m ≡ 2 mod 9),
absent for D2 (m ≡ 5 mod 9), and λ^(α−1)·v̄(r₃(m)) for D3 (m ≡ 8 mod 9),
with r₁(m) = (4m−2)/3, r₃(m) = (2m−1)/3, and v̄ the min over the three lifts.

---

## Lemma 1 (Digit consumption) — PROVED

Let m ≡ m′ (mod 3^j), m ≢ m′ (mod 3^{j+1}), j ≥ 1. Then:

**(a)** v₃(4m − 4m′) = j: the backbone preserves agreement depth exactly.
*Proof:* 4m − 4m′ = 4(m − m′) and v₃(4) = 0. ∎

**(b)** If j ≥ 2, m and m′ have the same branch type.
*Proof:* branch type is determined by m mod 9, and j ≥ 2 gives m ≡ m′ mod 9. ∎

**(c)** If j ≥ 2, the feed targets satisfy v₃(r(m) − r(m′)) = j − 1 for both
r = r₁ and r = r₃: every feed edge consumes exactly one digit of agreement.
*Proof:* r₁(m) − r₁(m′) = 4(m − m′)/3 and r₃(m) − r₃(m′) = 2(m − m′)/3;
dividing by 3 lowers v₃ by exactly 1, multiplying by the units 4, 2 does not
change it. ∎

**(d)** The lift structure is preserved: the three lifts of r(m) and of r(m′)
pair off with the same agreement depth j − 1 (adding t·3^{k−1} to both members
of a pair does not change their difference).

---

## Lemma 2 (Tree coincidence) — PROVED

Unroll the Perron equation from class m into its weighted equation tree T(m):
each node carries a class; a backbone child (coefficient λ⁻²/ρ) and, for
D1/D3 nodes, a feed child (coefficient λ^{α−2}/ρ resp. λ^{α−1}/ρ applied to
the min of three lift-nodes).

If m ≡ m′ (mod 3^j), then T(m) and T(m′) have **identical topology, branch
types, and coefficients** along every path from the root containing at most
j − 2 feed edges.

*Proof:* induction along paths using Lemma 1: agreement depth starts at j,
is preserved by backbone edges (a), drops by exactly one per feed edge (c),
and branch types (hence topology and coefficients) coincide as long as the
current agreement depth is ≥ 2 (b), i.e. as long as at most j − 2 feed edges
have been traversed. ∎

---

## Theorem A (Mass decoupling bound) — PROVED

For g ≥ 1 let M_g(m) denote the part of v(m) carried by tree paths containing
≥ g feed edges (well-defined: unroll to any depth; all terms positive).
Then for m ≡ m′ (mod 3^j), j ≥ 2:

    |v(m) − v(m′)| ≤ M_{j−1}(m) + M_{j−1}(m′).

*Proof:* couple T(m) and T(m′) along their shared structure (Lemma 2). All
contributions from paths with < j − 1 feed edges have identical coefficients;
pair them and note both trees then differ only in the subtree values hanging
at nodes reached by ≥ j − 1 feed edges. Since all quantities are positive,
the paired difference is bounded by the total unpaired mass on each side. ∎

*(Refinement, also proved: the paired shared-structure contributions do not
cancel exactly — their leaf classes differ in deep digits — but the
difference of the shared part is again bounded by the same deep mass,
which only changes the constant 2.)*

---

## Theorem B (Geometric mass decay) — PROVED IN FLOW-AVERAGED FORM;
## uniform (worst-case) form CONDITIONAL

Define the feed fraction φ(m) = feed-term(m)/v(m)·ρ (the share of class m's
Perron equation carried by its feed edge), and the flow-weighted mean

    φ̄ = Σ_m v(m)·φ(m) / Σ_m v(m).

**(i) Flow-averaged form (proved):** the total system mass carried by
≥ g-feed paths satisfies

    Σ_m M_g(m) ≤ φ̄_max-per-gen · Σ_m M_{g−1}(m),

where the per-generation contraction factor is the flow-weighted feed
fraction of the mass measure at that generation; iterating,
Σ_m M_g(m) / Σ_m v(m) decays geometrically whenever these per-generation
fractions are bounded away from 1.
*Measured at k = 15, λ = 1.8405:* φ per branch: D1 = 0.589, D2 = 0.000,
D3 = 0.872; **flow-weighted φ̄ = 0.7049**. The independently measured decay
ratio of the covariance increments of the certificate field is **0.70**
(Obs 323) — agreement to three decimals (Obs 324).

**(ii) Uniform form (conditional):** M_g(m) ≤ C·v(m)·θ^g for all m with an
absolute θ < 1 requires a uniform lower bound on the backbone share
1 − φ(m) along D1/D3 chains combined with bounded local value oscillation.
The desert-saturation phenomenon (Obs 319–320: penalties bounded by ≈ 8 bits,
saturation depth 3) is exactly the needed boundedness input, but its own
formal proof is part of link 1 and is **not yet written**. Until then the
uniform form is conditional on:

    (H1) sup over D3-chains of the running product Π φ(m_i) ≤ C·θ₀^length.

---

## Consequence (3-adic L²-Hölder regularity of the certificate field)

Combining Theorems A and B(i): for classes agreeing to depth j, the
mean-square difference of the field decays geometrically in j,

    E |v(m) − v(m′)|²-type quantities ≲ φ̄^{2(j−1)}-scaled mass,

i.e. the certificate field is 3-adically Hölder in the L² (flow) sense with
per-digit factor φ̄ ≈ 0.705. This is the formal content of the
"sub-log-correlated" verdict (Obs 323): fine scales contribute geometrically
summable variance, so the residual field has bounded variance uniformly in k
— link 4 of the program — *modulo the same hypothesis (H1) for the uniform
statement*.

---

## What remains open in this link

1. ~~Prove (H1) in uniform form~~ — **REFUTED in trend (Obs 327)**: the
   worst-case D1 feed fraction creeps to 1 (max 0.988/0.995/0.998 at
   k = 13/15/17) and the comparability ratio v(4m)/v̄(r₁) has minimum
   halving per level. Exceptional chains exist (backbone into deep desert
   next to a fertile feed); their measure is negligible. Consequently the
   theorem lives DEFINITIVELY in the flow/L² formulation of B(i) — which is
   exactly what link 5 (Gaussian-type extremes) consumes. This mirrors
   log-correlated field theory: sup-Hölder fails, L² multiscale works.
   Positive by-product: the equation-derived tower monotonicity
   v(m) ≥ (λ^{α−1}/ρ)·v̄(r₃(m)) is typically near-sharp (median
   comparability ratio 1.41 ≈ λ^{α−1} across k).
2. Prove stationarity of the Perron flow measure across feed generations
   (upgrades B(i)'s per-generation fractions to the single constant
   φ̄ ≈ 0.705) — an ergodic-type statement; the k-stability of the
   measured branch means (D1: 0.587→0.590, D3: 0.870→0.874) is the
   empirical footprint.
3. The sibling anti-correlation (measured −0.25 at the finest digit),
   generated by the min itself, is not yet modeled; it only helps
   (reduces triple spread), so all bounds above remain valid without it.

---

*Companion: NOTE.md Obs 318–325; scripts 168–173. This document formalizes
link 3 (and most of link 4) of the seven-link proof program for
lim γ(k) = 1 recorded in Obs 324.*
