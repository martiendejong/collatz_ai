# F16 — The Martingale Identity: E[3^k/2^(k+w)] = 1 exactly; excursion law P(>s) = 2^(−s)

**Status: PROVEN (identity) + CONFIRMED (tail slope −0.978 ≈ −1; extreme value 16.73 vs predicted 16.6).**
Figure: `figures/e16_extremes.png`

## The identity
Per macro-step the value multiplies by 3^k/2^(k+w). Under the (empirically exact, F12) geometric(½) law:

> E[(3/2)^k] = Σ (3/4)^k = **3**,  E[2^(−w)] = Σ 4^(−w) = **⅓**,  product = **1 exactly.**

**The Collatz value process is a martingale under the stochastic model** — log-scale drifts down
(−0.83 bits/macro-step) while linear-scale expectation is exactly conserved. The two facts coexist
because the multiplier distribution is heavy right-tailed.

## Consequence: the excursion law (optional stopping)
P(orbit ever exceeds 2^s × its start) ≤ 2^(−s), and ≈ 2^(−s) in fact.
Measured: tail slope **−0.978** bits per bit (predicted −1); max excursion among 100,000 orbits
**16.73 bits** (predicted log₂10⁵ = 16.6). No fat tail. (Note: an earlier naive prediction of −3.42
confused the single-step spike tail with the walk maximum; the correct exponent is the Lundberg
root θ*=1, which the martingale identity makes exact. Correction preserved for honesty.)

## Reading
- Explains the classical observation that record trajectory maxima scale like n² (excursion s ≈ log₂n).
- The martingale + optional stopping is the cleanest "proof" that divergence is impossible —
  under the i.i.d. model. Every constant is now exact; the sole unproven ingredient is that real
  orbits cannot deviate from the model forever (= E★, F15).

Related: [[10-alpha-star-theorem]], [[12-reload-independence]], [[15-borrow-chain-markov]]
