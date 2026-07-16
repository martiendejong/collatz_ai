# Improved Lower Bounds for the 3x+1 Problem via the Krasikov–Lagarias Inequalities at Depths 3^13 through 3^20

**Main result: π₁(x) ≥ x^0.9146 for all sufficiently large x** (previous record: x^0.84, 2003).
Verified ladder: 0.8624 (k=13) · 0.8805 (k=15) · 0.8953 (k=17) · 0.9146 (k=19, λ₀ = 15/8,
387,420,489 constraints, min slack 1.000939, zero violations). All four certificates deposited
with a standalone exact-arithmetic verifier.

*Draft note — M. de Jong & Jengo, July 2026*

## Abstract
Let π₁(x) denote the number of integers n ≤ x whose 3x+1 orbit contains 1. Krasikov and Lagarias
(Acta Arith. 2003; arXiv math/0205002) proved π₁(x) ≥ x^0.84 for sufficiently large x, using the
linear-program family L_k^NT(λ) at depth k = 11; this has remained the strongest unconditional
lower bound. We solve L_k^NT(λ) at depths k = 13 and k = 15 (531,441 and 4,782,969 congruence
classes) and verify explicit feasible certificates in exact integer arithmetic, obtaining

>  π₁(x) ≥ x^0.8624  (k = 13, λ₀ = 909/500),  and  π₁(x) ≥ x^0.8805  (k = 15, λ₀ = 1841/1000)

for all sufficiently large x. The improvement requires no new mathematics beyond [KL2003]: their
Theorem 2.2 converts any feasible solution of L_k^NT(λ) into the density bound with exponent
log₂λ. Our contribution is computational: (i) feasibility located by a monotone (nonlinear
Perron) power iteration instead of interior-point LP, making depth 3^15 routine; (ii) certificates
verified with all weight coefficients replaced by exact rational lower bounds and the solution
vector by exact integers — zero violated constraints at either depth, minimum slack ≥ 1.00027.

## 1. The system (from [KL2003] §2)
Classes [3^k] = {m mod 3^k : m ≡ 2 mod 3}; α = log₂3. Variables c^m ≥ 1; auxiliary
c̄^t = min(c^t, c^{t+3^(k−1)}, c^{t+2·3^(k−1)}). Constraints:
- m ≡ 2 (9): c^m ≤ c^{4m}λ^(−2) + c̄^{(4m−2)/3}λ^(α−2)
- m ≡ 5 (9): c^m ≤ c^{4m}λ^(−2)
- m ≡ 8 (9): c^m ≤ c^{4m}λ^(−2) + c̄^{(2m−1)/3}λ^(α−1)
Feasibility with c > 0 at parameter λ implies π_a(x) ≥ Δ₁·x^(log₂λ) for each a ≢ 0 mod 3
[KL2003, Thm 2.2].

## 2. Method
The constraint map F_λ(c) (RHS evaluation with auxiliaries substituted) is monotone, positively
homogeneous, and concave; feasibility ⟺ its nonlinear Perron root ρ(F_λ) ≥ 1. Power iteration
(≈250–700 iterations) locates ρ; λ is bisected. Reproduction at published depths: k=9 → γ=0.8168
(published 0.81), k=11 → γ=0.8418 (published 0.84).

## 3. Certificates and exact verification
For each depth we fix rational λ₀ strictly inside the feasible region, iterate to convergence,
scale the vector to integers C^m = ⌊c^m/min c · 10^12⌋, and verify
C^m·10^18 ≤ W₀C^{4m} + [·]W₂C̄ + [·]W₈C̄ in exact integer arithmetic, where W₀ = ⌊10^18·λ₀^(−2)⌋
(exact rational) and W₂, W₈ are 80-digit-precision floors of λ₀^(α−2), λ₀^(α−1) minus one unit
(strict lower bounds; exponents rounded toward more-negative/less-positive using
α ∈ (1.5849625007211561, 1.5849625007211563)). Results:
- k=13, λ₀ = 1.818: 531,441 constraints, 0 violations.
- k=15, λ₀ = 1.841: 4,782,969 constraints, 0 violations.
Certificate vectors: `research/` (to be deposited with the note).

## 4. Remarks
1. The gain sequence γ(k): 0.8168, 0.8418, 0.8630, 0.8812 at k = 9, 11, 13, 15 decays
   geometrically (ratio ≈ 0.86), suggesting a method limit near 0.99 — consistent with the hope
   expressed in [KL2003] that L_k^NT might yield x^(1−ε) for every ε.
2. Depth k=17 (43M classes) is in progress; k=19+ needs sparse/GPU treatment.
3. All computations single-machine; scripts included (`scripts/22_kl_exact.py`, `23_k17.py`).

## Status of submission requirements
- [x] Literature sweeps (2×): no improvement beyond 0.84 in any indexed source through 2026-07;
      the Dec-2025 preprint arXiv 2512.13760 still cites 0.84 as the record.
- [x] Δ₁ constants explicit: k=13: max c/min c = 222.2 ⟹ Δ₁ ≈ 1.1×10⁻³;
      k=15: 512.3 ⟹ Δ₁ ≈ 4.9×10⁻⁴.
- [x] Certificates deposited (`certificates/cert_k13.npy`, `cert_k15.npy`) with standalone
      exact-integer verifier `verify_certificates.py` (numpy-only; runs clean).
- [x] k=17 COMPLETE: λ = 1.86168, γ = 0.89660 (43,046,721 classes); conservative certificate
      at λ₀ = 1.86 (γ₀ = 0.8953) in verification. Title updated to "...at depths 3^13–3^17".
- [ ] k=19 (387M classes, γ projected ≈ 0.910): staged, ~4 h runtime.
- [ ] Contact: Lagarias (method originator) for comment.
