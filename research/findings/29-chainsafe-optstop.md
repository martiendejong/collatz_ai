# F29 — Three Semantic Bugs, One Honest Engine: γ ≥ 0.6527 (chain-safe, optimal stopping)

**Status: rigorous-form certificates (float bisection; rational certification pending).**
Script: inline round-13 (to be consolidated as `scripts/21_optstop.py` logic).

## The three bugs — each a mathematical lesson
1. **Refinement-SUM (tautology, γ→1):** summing u over trit-refinements is valid for CLASS COUNTS
   but not for a single leaf integer's subtree — it silently proves the heuristic, not a bound.
2. **Refinement-MIN per step (collapse, γ→0.5-floor):** worst case over unknown trits at every
   step gives an imaginary adversary compounding power no fixed integer has.
3. **Chain leak (round 12's 0.60–0.68 UNRIGOROUS):** harvesting info-exhausted leaves at full
   value assumes non-chain; an exhausted leaf could be ≡ 0 (mod 3) with a log-size subtree.
   Round 12's F28 numbers are hereby marked calibration-only, not certificates.

## The honest engine
Harvest only leaves with (i) residue info remaining (o < K), (ii) certified non-chain bottom trit
(r ≢ 0 mod 3), (iii) θ = 3^o/2^e < 1 — and use OPTIMAL STOPPING (any node may be declared a leaf;
max(stop, continue) DP), which is rigorous and recovers the exhaustion waste. All roots in one
backward DP — no sampling.

## Certified curve (worst root, full sweep)
| K | γ ≥ |
|---|---|
| 3 | 0.5000 |
| 4 | 0.5490 |
| 5 | 0.5837 |
| 6 | 0.6083 |
| 7 | 0.6309 |
| 8 | **0.6527** |

Monotone, ≈ +0.022/level. Naive extrapolation: K=12 → ~0.73. The remaining distance to
Krasikov–Lagarias 0.84 is their CLASS-COUNT difference-inequality system (π_a(x) functions at
mixed moduli, where the refinement-sum IS legitimate, LP-amortized) — the correct semantics
identified this round after the subtree-based sum proved tautological.

## Next (A1 continuation)
1. Reimplement with class-count semantics: π_a(x) ≥ π_{2^(-1)a}(x/2) + [parity] π-terms at
   mixed moduli with legit sum-consistency; LP over the inequality graph.
2. Rational/interval-arithmetic certification of the bisection output (theorem-grade).
3. Depth push K = 10–12.

Related: [[28-kl-engine-validated]] (superseded in rigor), PLAN.md Track A1
