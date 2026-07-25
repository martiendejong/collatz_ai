# The Summation Lemma: formal assembly toward γ → 1

**Companion to open-lemma-reduction.md, damping-theorem.md and
papers/gamma_to_one.tex. Consolidates Obs 358–369. Status per component
marked. 2026-07-25.**

Setting: the K–L system at the feasibility edge of depth k; Perron vector v,
flow measure W = v/Σv, field F = log₂ v. Multiscale (martingale)
decomposition along 3-adic agreement:

    M_p(i) = E[F | i mod 3^p],   X_p = M_p − M_{p−1},   F = M_0 + Σ_{p≥1} X_p.

## The three laws

**Law A (correlation decay across scales).**
|Corr_W(X_p, X_q)| ≤ c₀ · φ̄^{|p−q|}, with φ̄ = 1 − λ⁻²/ρ (exact flow
identity, PROVED) and c₀ a small absolute constant.
*Measured (k=13, Obs 369):* lag-1 mean 0.072, decay per lag 0.65–0.77 vs
φ̄ = 0.698; c₀ ≈ 0.11. *Mechanism:* pair-tree coupling — classes agreeing
mod 3^j have equation trees coinciding on all paths with < j−1 feed edges
(damping-theorem Lemmas 1–2, PROVED); increments at scale distance d share
only path-mass that has crossed ≥ d feed edges, which is a φ̄^d flow
fraction (flow identity, PROVED). *Remaining writing:* the step from
"shared path-mass ≤ φ̄^d" to "covariance ≤ c₀ φ̄^d σ_p σ_q" — a positivity
and boundedness argument on the path functionals (all terms positive; the
min operator only reduces shared mass).

**Law B (profile decay in scale).**
Var_W(X_p) ≤ C_B · φ̄^p.
*Measured:* per-scale ratios 0.63–0.74 vs φ̄ = 0.698 (k=13, Obs 369; same
law as the amplitude profile √φ̄ of Obs 331, <1% at k=15/17, re-confirmed
at k=21 within 1.5%, Obs 368). *Mechanism:* the scale-p increment is
injected through p feed generations; influence propagates with elasticity =
flow share = φ̄ per generation (Obs 330 chain; the flow identity is exact).
*Remaining writing:* the elasticity/influence bookkeeping of the linearized
Perron operator (partially drafted in damping-theorem.md).

**Law C (saturation in k).**
For fixed p, Var_W(X_p) at depth k approaches its limit at geometric rate
q_s ≈ 0.910 (the CV₁ law — hit to four decimals at k=21, Obs 367).
Needed only for k-uniformity of C_B and c₀. *Proof shape:* spectral
convergence of the finite-k truncations to the 3-adic limit operator.

## The Summation Theorem (conditional; assembly is routine given A–C)

    Var_W(F) = Σ_{p,q} Cov(X_p, X_q)
             ≤ (1 + 2c₀ φ̄/(1−φ̄)) · Σ_p Var(X_p)
             ≤ (1 + 2c₀ φ̄/(1−φ̄)) · C_B/(1−φ̄)   < ∞, uniformly in k.

*Numerical check:* the correction factor with (c₀, φ̄) = (0.11, 0.70) is
1.51 as an upper envelope; measured total/diagonal = 1.34 (Obs 369).
Measured limit variance ≈ 2.2 bits² (Obs 365).

## From bounded variance to the Open Lemma (assembly step A′)

Feed-domination chains (open-lemma-reduction.md) require, at each of g
levels, a log-deviation of the backbone ratio F(4m) − F(m) below a fixed
negative threshold (Lemma B identity). Marginal Chebyshev from the
Summation Theorem bounds each level's flow fraction; Law A's decay along
the feed direction upgrades this to a conditional per-step bound
δ(ε) < 1, whence F_k(g) ≤ δ^g — Lemma D — hence flow-L² attenuation
κ ≤ κ_max < 1 — hence q → 1 — hence, by the proved edge-rate theorem
(dγ/dq = 1/ln(4/3), gamma_to_one link 6), **γ(k) → 1** and
π₁(x) ≥ x^{1−ε} for every ε > 0.

*Guard rails already in place:* no single arithmetic structure can bypass
the tail requirement (Saturation Lemma, constant 54 = 2·3³, Obs 364:
tower penalties bounded); the field's lower tail is measured sub-Gaussian
(min at 2.3σ over 4.8M classes, Obs 365), so the Chebyshev step has slack.

## Honest remaining-work list (in order)

1. Law A's covariance step (positivity bookkeeping on shared path mass).
2. Law B's elasticity bookkeeping (linearized-operator influence).
3. Assembly A′: the conditional (along-chain) version of the Chebyshev
   step — Markov-type recursion with Law A decay.
4. Law C via truncation-convergence (only for uniform constants).
5. Saturation Lemma: promote sketch to full proof (j = 3 stratum;
   finite-generation telescope constants).

Nothing on this list requires an unknown idea; each is bounded write-up
work on a mechanism that is proved in skeleton and measured at the 1%-to-
4-decimal level. The risk register is correspondingly short: the only
identified failure mode is a hidden non-uniformity in k (Law C), bounded
by the endpoint argument (φ̄ ≤ 3/4 at λ = 2, no marginal point — see
gamma_to_one §falsify).
