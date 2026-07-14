# F12 — Reload Independence: the stochastic model is empirically exact, and that IS the wall

**Status: MEASURED (159,634 successive macro-step pairs along 4,000 orbits, 50–60 bit starts).**

## Results
- Correlation(k_i, k_{i+1}) = **−0.0026** (statistically zero at this n)
- Correlation(k_i, k_{i+2}) = +0.0051 (zero)
- Memorylessness: P(k′=1 | k) = .4992, .5012, .5078, .5076, .4905, .4793 for k = 1..6 — flat at ½
- χ² of reload law vs geometric(½): 3.18 (df 4, p ≈ 0.5) — indistinguishable from ideal

## Interpretation
Successive reloads along real orbits behave as i.i.d. geometric(½) draws to measurement precision.
Combined with F10's constants (drift −0.83, spike tail 2^(−3.42s)), the stochastic model predicts:
- P(orbit ever rises G bits above start) ≈ 2^(−cG) → no divergence, almost surely
- expected total stopping time ~ log n · const → Collatz behavior exactly as observed

**The paradox that defines the problem:** the reloads are fully deterministic (F1: computable from
(a mod 2^j, k mod 2^(j−2))), yet statistically indistinguishable from independent coin flips.
Every measurable trace of structure that could either prove or disprove long-run descent has
vanished into equidistribution. A proof must find structure that is invisible to correlation
statistics — or prove a form of equidistribution strong enough to cover every orbit
(which is precisely the ÷2-in-ternary equidistribution frontier, F2).

Related: [[10-alpha-star-theorem]], [[09-verdicts-and-open-core]]
