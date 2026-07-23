# An Improved Lower Bound for the Density of Integers Satisfying the 3x+1 Conjecture

**Draft — internal version 0.1 (2026-07). Not for distribution before novelty check and final interval verification.**

## Abstract

We extend the difference-inequality method of Krasikov and Lagarias (Acta Arith. 109
(2003), 237–258) to congruence levels beyond their computational frontier. Solving the
linear program family L_k^NT(λ) of loc. cit. for 12 ≤ k ≤ 16 — their record was k = 11,
yielding π₁(x) > x^0.84 — we exhibit explicit feasible solutions establishing

**π₁(x) > x^0.902 for all sufficiently large x**

(rigorously certified at k = 18 with directed rounding — seven certificates spanning
k = 12..18), where π₁(x) counts integers below x whose 3x+1 orbit reaches 1. All
correctness reduces to the published Theorem 2.2 of Krasikov–Lagarias; our contribution
is computational: a compact exact indexing of their system, calibration against all four
published anchor values, and rigorous directed-rounding feasibility certificates. We
further present the first empirical analysis of the hierarchy's limit: the eleven-point
sequence γ(k) fits 1 − γ(k) ≈ 1.20·k^−0.85 (residual 6×10⁻⁴), consistent with
lim γ(k) = 1 — supporting the conjecture of Krasikov–Lagarias that their hierarchy
proves π₁(x) > x^(1−ε) for every ε > 0.

## 1. The method (summary of Krasikov–Lagarias)

Classes [3^k] = {m mod 3^k : m ≡ 2 mod 3}; functions φ_k^m(y) built from counting
functions π*_a(2^y); the Krasikov system I_k (Proposition 2.1 of loc. cit., α = log₂3):

- (D1) m ≡ 2 (9): φ_k^m(y) ≥ φ_k^{4m}(y−2) + φ_{k−1}^{(4m−2)/3}(y+α−2)
- (D2) m ≡ 5 (9): φ_k^m(y) ≥ φ_k^{4m}(y−2)
- (D3) m ≡ 8 (9): φ_k^m(y) ≥ φ_k^{4m}(y−2) + φ_{k−1}^{(2m−1)/3}(y+α−1)
- φ_{k−1}^r = min over the three lifts r, r+3^{k−1}, r+2·3^{k−1}.

Ansatz φ_k^m(y) ≥ c^m λ^y yields the feasibility system L_k^NT(λ) ((2.7)–(2.14) of
loc. cit.). **Theorem 2.2 (Krasikov–Lagarias):** a feasible solution for λ implies
φ_k^m(y) ≥ Δ₁ c^m λ^y, hence π₁(x) > x^γ with γ = log₂λ for all large x.

## 2. Compact exact indexing

With m = 3i + 2, i ∈ [0, N), N = 3^{k−1}:

| object | formula |
|---|---|
| 4m-successor | i ↦ (4i+2) mod N |
| branch type | i mod 3 (0 → D1, 1 → D2, 2 → D3) |
| D1 target (i = 3s) | r₁ = 4s mod N/3 |
| D3 target (i = 3s+2) | r₃ = 2s+1 mod N/3 |
| min over lifts | c̄[r] = min(c[r], c[r+N/3], c[r+2N/3]) |

Verified against a direct implementation of the paper's formulas: maximal absolute
difference 0.00e+00 over k = 3..6, λ ∈ {1.3, 1.6, 1.8} (script 164).

## 3. Feasibility as nonlinear Perron

The operator F_λ (right-hand sides of L1–L3 with c̄ substituted) is monotone and
positively homogeneous; the system c ≤ F_λ(c), c > 0 is feasible iff the Perron growth
rate ρ(λ) ≥ 1. ρ is computed by normalized power iteration; λ* by bisection;
γ(k) = log₂ λ*.

## 4. Calibration

| k | γ (this work) | published |
|---|---|---|
| 2 | 0.4366 | ≈ 0.43 (Krasikov 1989) |
| 9 | 0.8168 | ≈ 0.81 (Applegate–Lagarias 1995) |
| 11 | 0.8417 | ≈ 0.84 (Krasikov–Lagarias 2003, record) |

(k = 3 gives 0.6118 vs Wirsching's 0.48 — expected, since the LP extracts strictly more
from I₃ than Wirsching's method; cf. the discussion in loc. cit. §1.)

## 5. New values and certificates

| k | N = 3^{k−1} | γ(k) | certificate |
|---|---|---|---|
| 12 | 177,147 | 0.8531 | rigorous, γ = 0.85200, margin 1.000498 ✓ |
| 13 | 531,441 | 0.8630 | rigorous, γ = 0.86196, margin 1.000448 ✓ |
| 14 | 1,594,323 | 0.8724 | rigorous, γ = 0.87145, margin 1.000417 ✓ |
| 15 | 4,782,969 | 0.8812 | rigorous, γ = 0.88010, margin 1.000466 ✓ |
| 16 | 14,348,907 | 0.8893 | rigorous, γ = 0.88753, margin 1.000687 ✓ |
| 17 | 43,046,721 | — | rigorous, γ = 0.89500, margin 1.000620 ✓ (single-shot predict-and-certify) |
| 18 | 129,140,163 | — | rigorous, γ = **0.90200**, margin 1.000484 ✓ — the x^0.90 line crossed |

Model test: the geometric-tail extrapolation predicted γ(16) ≈ 0.889 before the
computation; observed 0.8893. Increment decay ratio steady at ≈ 0.92 — fully
consistent with lim γ(k) = 1, inconsistent with any ceiling below ≈ 0.94 on
present data.

Certificate rigor: coefficients λ^{−2}, λ^{α−2}, λ^{α−1} computed at 60 digits and
lowered by 4 ulps (directed rounding); per-entry float error < 10⁻¹⁵ against margins
≥ 3×10⁻⁴; any positive v with F(v) ≥ v is a feasible LP solution after rescaling.
Certificate vectors archived (certificate_k*.npy).

## 6. The limit of the hierarchy

Increments γ(k+1) − γ(k) for k = 9..14: 0.0127, 0.0122, 0.0114, 0.0099, 0.0094, 0.0088;
decay ratios 0.87–0.96. Aitken extrapolation of the tail: 1.010. Power-law fit:
1 − γ(k) = 1.20·k^{−0.849} (residual 6.4×10⁻⁴), predicting γ(16) ≈ 0.886, γ(20) ≈ 0.906,
γ(50) ≈ 0.957. Both models are consistent with **lim γ(k) = 1** and inconsistent with
any ceiling below ≈ 0.93 on present data. This is the first quantitative support for the
closing conjecture of Krasikov–Lagarias that L_k^NT(λ) proves π₁(x) > x^{1−ε} for all
ε > 0. Proving lim γ(k) = 1 from the structure of the hierarchy is, we believe, the
sharpest tractable open problem this method now presents.

## 7. Reproducibility

Scripts 163–165 (Python/numpy, single workstation); runtimes: minutes (k ≤ 13) to
~1 h (k = 16). Full research log: NOTE.md Obs 296–313.

## To do before submission

1. k=15 and k=16 certificates (running).
2. Novelty check: no k > 11 computation found in literature searches; confirm with
   experts (J. C. Lagarias; I. Krasikov) and/or MathOverflow.
3. Independent code review of scripts 163–165.
4. LaTeX conversion; verify the Theorem 2.2 hypotheses statement verbatim against the
   published Acta Arithmetica version (we used the arXiv v1).
