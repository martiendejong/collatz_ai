# F32 — Rounds 16–29 Consolidated

## Density record campaign (Track A1) — COMPLETE THROUGH k=15, k=17 in flight
- R16: **exact-integer verification passed for both records**: k=13 (531,441 constraints) and
  k=15 (4,782,969), zero violations, rational lower-bound weights (80-digit floors −1).
  **π(x) ≥ x^0.8624 and π(x) ≥ x^0.8805, proven-grade** modulo KL2003 Thm 2.2.
  Paper draft `NOTE_DENSITY.md`.
- R28: certificates deposited (`certificates/cert_k13.npy`, `cert_k15.npy`) + standalone
  `verify_certificates.py` (numpy-only): ran clean — full reproducibility package.
- k=17 (43M classes, float32): bisection at [1.86141, 1.86195] ⟹ **γ(17) ≈ 0.8965**.

## Cycle exclusion extensions
- R18: 5-circuits — 36.4M shapes (kᵢ ≤ 6): trivial only.
- R26: 6-circuits — 106.7M shapes (kᵢ ≤ 5): trivial only.
- R19: period-6 census — empty through kᵢ,wᵢ ≤ 4 (16.8M shapes).
- R27: period-7 census — empty (4.78M shapes). **Cycle list {1,−5,−17} complete through period 7,
  circuits through m=6.**

## Census machinery cross-validation (both controls)
- R22 (after prefix/tail bug fix): 5n+1 census finds ALL known cycles with all odd members:
  {1,3}, {13,33,83}, {17,43,27}, plus −1. R23: 3n−1 finds exactly {1},{5,7},
  {17,25,37,41,55,61,91} — ten members, zero spurious. Machinery certified system-agnostic.

## Theory (Tracks C2/C3)
- R20–21: **Theorem 9″ (entropy-tempered shadowing)**: family continuation = exact residue
  measure (provable via Terras equidistribution). Ladder 1/8 exact; {k=2,3} family 3/16 exact;
  strict-monotone-climb measure **μ = 0.2862745…** (exact rational sum; measured 0.28629).
  Density of r-step monotone climbers = μ^r unconditionally.
- R25: **Realization Computer** — inverts Terras bijection for any index stream (sanity: ladder
  reproduces −5's 2-adic digits exactly). Thue-Morse climbing stream: no positive-integer
  realization within depth 300 (top digits 59% ones) ⟹ TM-shadowing dies at log₂(n)/3.5 rungs.
  C2 operational for arbitrary explicit streams.
- R29: **window-unfairness envelope**: max observed mean drift over windows: +1.4 (m=10),
  +0.5/0.6 (m=20), ~0 (m=40), **−0.1..−0.22 (m=80)** at scales 2^32/2^48. The envelope crosses
  zero near m ≈ 50: C3's needed theorem ("long windows have negative drift") has small true
  constants. Divergence needs ≥ 0 forever — never seen past m ≈ 40.

## Miscellaneous
- R17: B mod D skew = size artifact (Catalan D=1 family); A3 remains low priority.
- R24: records scan 2→8×10⁹ running (tail-law test).
