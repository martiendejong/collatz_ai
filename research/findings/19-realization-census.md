# F19 — The Realization Census: the integer-realized periodic streams are EXACTLY the known cycles

**Status: PROVEN (exhaustive in range) — the strongest structural result of the program.**

## Setup
Every periodic index stream (k₁,w₁)…(k_m,w_m) has a unique 2-adic realization
a₁ = B(shape)/(2^(K+W) − 3^K) (F13/F14 algebra, sign free). It is a genuine cycle of the Collatz map
on ℤ iff a₁ is an odd integer. Census over all shapes:

## Period 1 (all k, w ≤ 40)
| shape | a | n | identity |
|---|---|---|---|
| (1,1) | 1 | **1** | the trivial cycle 1→4→2→1 |
| (2,1) | −1 | **−5** | the negative cycle −5→−14→−7→−20→−10 (= infinite 9/8 ladder, F17) |

**No others.** 1,600 shapes; two integer realizations.

## Period 2, non-repeat (all kᵢ, wᵢ ≤ 18)
| shape | a | n | identity |
|---|---|---|---|
| (4,1),(3,3) | −1 | **−17** | the long negative cycle (18 steps) |
| (3,3),(4,1) | −5 | −41 | same cycle, second phase |

**No others.** 104,976 shapes; one cycle (in its two rotations).

## Circuits m ≤ 4, positive side (F13/F14/E18)
m=2 (shapes ≤ 60), m=3 (≤ 24), m=4 (≤ 10): only the trivial cycle. 32.5M shapes at m=4 alone.

## Why this matters
The census machinery **finds every cycle that exists in its range** — it independently rediscovered
1, −5, and −17 with no prior knowledge. The emptiness of the positive side is therefore a genuine
structural fact within the searched region, not a blind spot. (The degenerate −1→−2→−1 cycle sits
at the coordinate boundary n+1=0 and is outside the (a,k) chart.)

**The Collatz cycle problem, exactly restated:** the map shape ↦ B(shape)/(2^(K+W) − 3^K)
hits an odd POSITIVE integer only at the trivial shape. Known hits: +1, −1, −5 (as a₁ values
1, −1, −5 across phases). The conjecture's cycle half = "no further positive hits at any period."

Related: [[13-two-circuit-exclusion]], [[14-three-circuit-exclusion]], [[17-ladders-and-minus5]]
