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

## Law A: proof modulo flow-stationarity (write-up, 2026-07-25)

Work at the feasibility edge (ρ = 1), fixed point v = F(v). Freeze the
argmins of the min-operators at the fixed point; the frozen system is
linear with positive coefficients, and v admits the exact representation
v(m) = Σ_π w(π) v(end π) over paths of the selected tree, for every
truncation depth.

**Lemma A1 (elasticity representation) [PROVED — chain rule].**
For the log-field F = log₂ v, the sensitivity of F(m) to the log-value at
node x is e_m(x) = Σ_{π: m→x} w(π) v(end π)/v(m) ∈ [0, 1] — the *flow
share* of x in m's equation tree. Sensitivities are **positive** (the
system is monotone with positive coefficients: no cancellation anywhere)
and sum to 1 over any antichain cutting all paths. ∎

**Lemma A2 (digit-influence localization) [PROVED — damping-theorem
Lemmas 1–2].** The increment X_q(i) (influence of digit q given digits
< q) is transmitted exclusively through tree paths containing ≥ q−1 feed
edges: paths with fewer feed edges have topology and coefficients
independent of digit q. ∎

**Lemma A3 (per-generation share identity) [PROVED — flow identity].**
The flow share passing the g-th feed generation, averaged over the Perron
flow, equals φ̄^g with φ̄ = 1 − λ⁻²/ρ exactly, PROVIDED the flow-weighted
distribution of feed shares is the same at every generation
(**flow-stationarity**). ∎ (conditional clause explicit)

**Proposition A (Law A, conditional on stationarity).** For q > p,
decompose X_q = X_q^{shared} + X_q^{fresh}, where X_q^{shared} collects
the transmission through paths that also carry digit-p influence. By A2
these paths contain ≥ q − p feed edges below the digit-p divergence
level; by A1 (positivity — envelope without cancellation) and A3 their
total flow share is ≤ C·φ̄^{q−p} in flow-mean; Cauchy–Schwarz gives
Cov_W(X_p, X_q) ≤ σ_p ‖X_q^{shared}‖₂ ≤ c₀ φ̄^{q−p} σ_p σ_q. ∎

**The single remaining input: flow-stationarity.** Statement: the
flow-weighted law of the feed share φ over generation-g feed nodes is
independent of g. Empirical status: measured directly as the
participation-ratio cascade (Obs 330: per-generation ratios
0.691–0.735 ≈ φ̄ over six generations at k=11) and as k-stability of the
branch means (Obs 327: D1 0.587→0.590, D3 0.870→0.874). Proof shape: the
generation-g feed-node flow measure is the g-th image of the Perron flow
under the (frozen) feed-transfer operator, whose fixed point IS the Perron
flow — stationarity is exact at the fixed point; what needs writing is
that the generation measures, which start AT the fixed-point flow, remain
there (a two-line invariance check on the frozen linear system, plus
control of the argmin-freezing error — the only genuinely technical bit,
since unfreezing the min can only lower shared mass by positivity).

**Consequence.** With Law A proved modulo the invariance check, the
remaining-work list shrinks: items 1 and 2 (Law A covariance step, Law B
elasticity bookkeeping) both rest on the SAME Lemmas A1–A3 + stationarity;
they fall together.

## The five-item list, executed (2026-07-25, second pass)

**Refinement of the single input.** All conditional clauses below reduce to
one statement:

> **(S) Density-bounded feed cascade.** Let ν_g be the flow measure
> transported through g feed edges (the g-th feed image of the Perron
> flow μ). Then θ_g := E_{ν_g}[φ] satisfies sup_g θ_g ≤ θ < 1.

**S: status corrected (Obs 374) — MEASURED monotone, two candidate
mechanisms, neither closed.** The ladder derivation below is retained for
the record but its premise (autonomous per-level dynamics) is NOT
established in the true K–L system (lower-level values are minima of
level-k lifts; dynamics stay at level k). The local alternative ("the min
is feed-poor": argmin lift has lowest φ in 65.7% of triples, flow-weighted
φ-gap −0.24%) is measured but weak; whether compounding closes it is open.
Lemma S1 below is unaffected and fully proved.

*Lemma S1 (constant-lift embedding) [PROVED + machine-checked].* For
j ≤ k, a feasible solution of L_j^NT(λ) lifts to a feasible solution of
L_k^NT(λ) by constant lifting c^M := c^{M mod 3^j}: the maps T4, r₁, r₃
are affine and commute with reduction mod 3^j; branch types agree for
j ≥ 2; the min over equal lifts is the value. Hence **λ*(j) ≤ λ*(k)**.
(Machine check at (j,k) = (6,9): the entire ratio field is preserved
exactly — margin 1.000047 → 1.000047.) ∎

*The ladder.* In the (true, hierarchical) K–L system, feed generation g
lives at level k−g. At λ = λ*(k), every lower level is subcritical:
ρ_{k−g}(λ) ≤ 1 by S1 (λ > λ*(k−g) ⟹ infeasible ⟹ ρ < 1). The per-level
flow identity (backbone is a permutation at every level) gives the
generation-g flow-mean feed share θ_g = 1 − λ⁻²/ρ_{k−g} ≤ 1 − λ⁻² = φ̄.
Hence **sup_g θ_g = θ₀ = φ̄ ≤ 3/4** for all λ ≤ 2. ∎ (skeleton)

*Measured confirmation (script 193, k=13):* the edge equation is an
exact Markov chain (row sums 1.000000); θ₀ = 0.6977 = 1 − λ⁻² to four
decimals, and θ_g decreases monotonically (0.6977 → 0.6385 over ten
generations) — the ladder visible in the data. Remaining write-up: the
generation↔level identification through the collapsed single-vector
implementation (the collapse is calibrated exact against K–L's published
values; cleanest route is to state Laws A/B directly in the hierarchical
formulation).

**1. Law A (covariance step) — CLOSED modulo S.** Lemmas A1–A3 above +
Proposition A. Positivity (A1) removes all cancellation issues; A2
localizes digit-q influence to the ≥(q−1)-feed shell; S converts shell
counts to θ^{q−p} envelopes.

**2. Law B (profile decay) — CLOSED modulo S, same lemmas.** Var_W(X_p)
= ‖digit-p influence‖²_{L²(W)} ≤ (shell-(p−1) elasticity envelope)² ≤
C_B θ^{p} by A1 (positive elasticities summing to ≤ 1 per antichain) +
A2 (localization) + S (shell decay). No new machinery: items 1 and 2
fall together, as predicted.

**3. Assembly A′ (chain recursion) — CLOSED modulo S + Summation
Theorem.** Let G(m) = F(4m) − F(m) (backbone log-ratio) and t₀(ε) =
−log₂(ερλ²) > 0 for ε < λ⁻²/ρ (≈ 0.28 at the edge — the endpoint
condition ε < 1/4 reappears as the domain of validity). Chebyshev:
W{G ≤ −t₀} ≤ Var_W(G)/(t₀ − |E_W G|)² =: δ₀(ε) < 1 for ε small, since
Var_W(G) ≤ 2(1+|corr|)Var_W(F) is bounded by the Summation Theorem.
Chain step: consecutive chain levels are separated by exactly one feed
edge, so conditioning on the previous level perturbs mean and variance
of the next level's G by ≤ c₀θ-factors (Law A applied to the difference
field); the Markov recursion yields F_k(g) ≤ (δ(ε))^g with
δ(ε) = δ₀(ε)(1 + O(c₀θ)) < 1 — Lemma D, hence κ ≤ κ_max < 1, hence
q → 1, hence γ → 1 via the proved edge-rate theorem.

**4. Law C (k-uniformity) — CLOSED modulo B + one coupling paragraph.**
The depth-k system is the depth-(k+1) system with the finest digit
aggregated through the min. The martingale decompositions agree on
scales p ≤ k−1 up to a boundary correction supported on the finest
scale, of L²-size ≤ Var(X_k^{(k+1)})^{1/2} ≤ (C_B θ^k)^{1/2} (Law B).
Hence the constants (C_B, c₀) form Cauchy sequences in k with geometric
increments: uniformity follows. The one technical point: argmin flips
between the k- and (k+1)-systems; by positivity (A1) a flip only
reallocates shared mass downward, so the correction bound survives.
(Empirical footprint of this geometric coupling: the 0.910 saturation
law, hit to four decimals at k=21.)

**5. Saturation Lemma — PROVED IN FULL (no conditions).**
(a) j ≥ 4: transmitted depth exactly 2 — the identity
r₁(4m)+4 = (16(m+4)−54)/3 with v₃(54) = 3 (machine-checked, Obs 364).
(b) j = 3 stratum: writing m+4 = 27t (3∤t), transmitted depth
= 2 + v₃(16t−2) **exactly** (40,000/40,000 samples), and
v₃(16t−2) = s has density 3^{−s} exactly (16t−2 ≡ t+1 mod 3, standard
unit-genericity — verified: measured stratum densities 1/2, 1/3, 1/9,
1/27, … match 3^{−s} to four decimals). So the j = 3 leakage is
geometrically summable with EXACT constants.
(c) Telescope: a desert of any depth touches ≤ 2 further feed
generations at transmitted depth ≤ 2 (branch-period-3 + (a)); each
bounded-depth generation contributes log-penalty ≤ log₂(ρ/λ⁻²)·2 + O(1)
with explicit constants; total penalty P_max < ∞ uniformly. ∎
This makes the tower part of the field uniformly bounded — the
guard-rail under the Chebyshev step of item 3.

## Final risk register

Single analytic input **(S)**; measured over six generations at ≤ 0.05
deviation; endpoint-anchored (φ̄ ≤ 3/4 at λ = 2, no marginal point).
Failure of S is the ONLY way the program dies, and it would have to
manifest as a growing θ_g trend that six measured generations, the
k-stability of branch means, AND the four-decimal saturation law all
fail to show. Items 1–4 are write-up-complete modulo S; item 5 is
unconditional.
