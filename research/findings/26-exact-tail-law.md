# F26 — The Exact Tail Law: P(σ ≥ t) = C · t^(−3/2) · 2^(−t/20)

**Status: DERIVED (Cramér, under the coin model) + CONFIRMED flat over five orders of magnitude.**

## Derivation
A T-step contributes X = log₂3 − 1 (odd, prob ½) or −1 (even, prob ½). Survival (σ ≥ t) = the walk
stays ≥ 0 for t steps. The Cramér exponent is min_θ E[2^(θX)] = ½(2^(0.585θ) + 2^(−θ)), minimized
at θ = 0.488 with value **0.96594 = 2^(−0.0500)** — the needle KL rate, recovered analytically.
Random-walk theory adds the ballot prefactor t^(−3/2):

> **P(σ ≥ t) = C · t^(−3/2) · 2^(−t/20), C ≈ 15.  "One bit per twenty steps."**

## Confirmation (1.5M orbits, 30–34 bits)
Residual log₂P + t/20 + 1.5·log₂t: **3.35, 3.73, 3.83, 4.03, 4.09, 3.96, 4.41, 3.93** for
t = 20…250 — flat within ±0.5 bits while P falls from 5.7e−2 to 6.7e−7.

## The rate-bend of F25, explained exactly
Apparent (window) slopes: [20,140] → 0.05 + 1.5·Δlog₂t/Δt = 0.085 (measured 0.080 ✓);
t ≈ 350 → 0.056–0.06 (records implied 0.064 ✓). There is no scale-dependent physics —
just the polynomial prefactor riding the pure exponential 2^(−t/20). F25's falsification is
fully resolved: the linear records fit failed because it absorbed the prefactor into the rate.

## Records check
Expected max solves log₂(M·C) = 1.5·log₂t + t/20. Scan complete to 2×10⁹; true record chain tail:
63728127:376 (−0.7 bits vs EV) · 217740015:395 (0.0) · 1200991791:398 (−2.3) ·
**1827397567:433 (+2.1, a 25%-probability fluctuation)**.
The whole record chain — 27 through 1.8×10⁹ — sits on one curve with a single fitted constant
(the rate 1/20 and exponent −3/2 are parameter-free theory).

## Status note
The law is exact under the coin model; its per-orbit validity is H(ε)-conditional like everything
in the bulk. As an instrument it converts any future record scan into a sharp hypothesis test of
the model at depth t — the deepest available probes of E★.

Related: [[25-records-falsification]], [[07-never-drop-needle]], [[24-theorem-c-conditional]]
