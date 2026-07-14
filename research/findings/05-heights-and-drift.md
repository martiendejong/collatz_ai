# F5 — Height Functions: uniform drift −0.83, no monotone weighting exists

**Status: MEASURED (120,000 macro-steps, 30–46 bit inputs; blocks 4,000 × 40–60 bit).**

## The α-sweep: H_α(a,k) = log₂(a) + α·k over one macro-step

| α | mean ΔH | P(increase) | max increase |
|---|---------|-------------|--------------|
| 0.0 | −0.827 | .327 | 25.1 |
| 0.585 | −0.829 | .287 | 15.8 |
| 0.8 | −0.830 | .275 | 12.3 |
| **1.2** | −0.832 | **.242** | **5.9** |
| 1.585 | −0.834 | .283 | 8.4 |
| 2.0 | −0.836 | .285 | 14.6 |

- Mean drift is **α-independent** (≈ −0.83 bits/macro-step) because E[Δk] ≈ 0 (k is geometric-stationary). The drift lives entirely in the family parameter.
- No α is monotone. Best excursion control near **α ≈ 1.2** (increase probability .24, max +5.9) — candidate weighting for future work; not obviously log₂3 or log₂3−1.

## Blocks of B macro-steps (α = 0.585)
| B | P(H increases over block) | max increase |
|---|--------------------------|--------------|
| 1 | .274 | 9.8 |
| 5 | .134 | 12.2 |
| 10 | .054 | 9.9 |
| 20 | **.0088** | **11.2** |

**The key asymmetry: probability collapses with block length, magnitude does not.** Rare-but-large excursions are exactly the "almost all vs all" gap in numerical form. Any proof must bound magnitude, not just probability.

## Control
5n+1 with the same H: mean **+0.314**, P(inc) = .715. Opposite sign — the height battery passes the litmus test; every drift statement genuinely uses log₂3 < 2.

Related: [[06-ternary-streaks]], [[09-verdicts-and-open-core]]
