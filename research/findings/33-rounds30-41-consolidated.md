# F33 — Rounds 30–41 Consolidated

## Density campaign
- R31/33: Δ₁ constants explicit (k=13: ≈1.1e−3; k=15: ≈4.9e−4) — bounds fully explicit.
- R32: literature sweep #2 — nothing beyond 0.84 anywhere through 2026-07. Results novel.
- R34/35: RAM audit (31.3 GB) — k=19 GO (needs ~10 GB); script staged; γ(19) projected 0.9100.
- R39: γ-ceiling: gain ratios stable at 0.854; two extrapolations (geometric, Shanks) agree:
  **method ceiling γ∞ ≈ 0.985–0.987**. Reaches x^0.97–0.99, not x^1 — the final density gap
  needs the k→∞ structure theory K–L hoped for.
- k=17 bisection: [1.86168, 1.86182] ⟹ γ(17) ≈ 0.8967 (final line pending).

## Model tests (all passed)
- R40: records scan COMPLETE to 8×10⁹: final record 2788008987:σ=447, on the tail-law curve.
  One record expected in range per the law; exactly one found.
- R38: excursion law at 2^60: slope −1.018 (theory −1) — martingale scale-invariant.
- R41: borrow bias dev = 0.8/bitlength stable 24→224 bits — boundary layer only, no bulk
  structure at any tested depth.
- R36: Chang's one-bit statistic unbiased to 1.6e−4 — all known reductions rest on
  empirically-perfect hypotheses; the wall is shared.

## Theory
- R37: Theorem 10 (entropy-tempered shadowing, μ = 0.2862745 corollary) added to NOTE.md.
- R29 (in F32): window-unfairness envelope crosses zero at m ≈ 50.

## Program state at R42
Ten theorems + two propositions (NOTE.md); density record 0.84 → 0.8805 verified (0.8966, 0.910
in flight); cycles empty through period 7 / circuits m ≤ 6; census machinery cross-validated on
both control systems; tail law tested to depth 447 and scale 8×10⁹; reproducibility package
deposited. The conjecture's open core (H(ε)/E★ + Baker-throttled large-period cycles) unchanged —
every new measurement lands exactly where the theory says it must.
