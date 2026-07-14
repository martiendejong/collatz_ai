# F10 — The α★ Theorem: excursion-optimal height weight is (1+log₂3)/2

**Status: PROVEN (elementary) + confirmed on 200,000 macro-steps.** Figure: `figures/e10_alphastar.png`

## Statement
For the height family H_α(a,k) = log₂a + α·k, the macro-step change is (exactly, up to O(1/n)):

> ΔH = k(log₂3 − α) − w + k′(α − 1)

Large spikes come from either a large current index k (slope log₂3 − α) or a large reload k′ (slope α − 1).
The worst-case spike slope is minimized where the two balance:

> **α★ = (1 + log₂3)/2 = 1.29248…, spike slope (log₂3 − 1)/2 = 0.29248…**

## Confirmation (fine sweep, 200k samples, 40–58 bit inputs)
Max increase: 10.7 (α=1.0) → declining linearly → **5.9 at α=1.3** → rising linearly → 10.9 (α=1.65).
Clean V with vertex at α★. Worst observed spikes match the slope law: ΔH=6.02 at k=20 (0.2925·20=5.85);
ΔH=5.73 at k′=18 (0.2925·18=5.26 + cascade term).

## Why it matters
At α★, H has: drift **−0.83 bits/macro-step**, spike tail **P(spike > s) ≈ 2^(−s/0.2925) = 2^(−3.42·s)**
(both index tails geometric(½)). These are the sharpest Foster–Lyapunov constants of the program:
a biased walk with drift −0.83 and exponential-3.42 increments "cannot" diverge — except the increments
along one orbit are deterministic, not independent. The independence gap (F12) is all that separates
these constants from a proof of no-divergence.

Related: [[05-heights-and-drift]], [[12-reload-independence]]
