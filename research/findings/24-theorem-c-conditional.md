# F24 — Theorem C: the Collatz Conjecture under Quarter-Fair Mixing

**Status: PROVEN conditional theorem. This names the exact assumption the conjecture rests on.**

## Setup
Macro-step i of an orbit has index kᵢ and cascade wᵢ; the exact value change is
log₂(mᵢ₊₁/mᵢ) = kᵢ·log₂3 − kᵢ − wᵢ + ηᵢ with |ηᵢ| < 2^(−69) whenever mᵢ > 2^71.
At the residue-class level it is PROVEN (units act bijectively mod 2^j) that
P(k ≥ t) = 2^(1−t) and P(w ≥ t) = 2^(1−t) exactly over full classes.

## Hypothesis H(ε) — "quarter-fair graded mixing"
Along every orbit that eventually stays above 2^71 (if any exists), the asymptotic frequencies obey
1. **freq(kᵢ ≥ t) ≤ (1+ε)·2^(1−t) for every t ≥ 1** (graded reload fairness), and
2. **liminf mean of min(wᵢ, 8) ≥ (1−ε)·(2 − 2^(−7))** (truncated cascade fairness).

## Theorem C
**If H(ε) holds for some ε < ¼, the Collatz conjecture is true.**

*Proof.*
(1) Mean k: Σₜ freq(k ≥ t) ≤ 2(1+ε). Mean w ≥ mean(w∧8) ≥ (1−ε)(2−2^(−7)).
Mean log₂-change per macro-step ≤ (log₂3−1)·2(1+ε) − (1−ε)(2−2^(−7)) + 2^(−69)
= **−0.8223 + 3.1621·ε** < 0 for ε < 0.26 (at ε = ¼: −0.0317 < 0).
(2) *No divergence:* a divergent orbit must eventually stay above 2^71 (each return below enters
the verified region and converges — contradiction), so from some point its partial log-sums drift
to −∞ at rate ≥ 0.03 bits/macro-step, forcing the value below 2^71 — contradiction.
(3) *No large cycle:* around a cycle the mean log₂-change is exactly 0, contradicting the bound.
Cycles with min element ≤ 2^71 are excluded by verification except {1,4,2}. ∎

## What this achieves and what it does not
- The divergence half AND the cycle half both reduce to ONE hypothesis about visiting frequencies
  of residue classes — a mixing statement, not a drift statement (drift is derived, not assumed).
- The tolerance is enormous: fairness up to a factor 1.25 suffices.
  **Measured reality: ε ≈ 0.0002 (cascades), ≈ 0.01 (reloads at statistically reliable depths).**
  The safety margin between measured and required is factor ~25–1000. The margin between measured
  and PROVEN is the conjecture.
- H(ε) is exactly the per-orbit equidistribution of the ×3-plus-shift action on 2-adic residues —
  the same genre as the open Furstenberg ×2×3 rigidity problems. Theorem C does not evade the wall;
  it gives the wall its precise name and its price: **19–26% fairness, forever, for every orbit.**

Related: [[15-borrow-chain-markov]], [[16-martingale-identity]], [[23-synthesis-consistency-web]]
