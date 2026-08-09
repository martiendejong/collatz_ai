# 364: C4 from the 2026-08-09 analysis (the never-written "script 175"):
# which rate governs climb sojourns? Obs 180 gives the naive Cramer rate
# I(0) = 0.2113 bits/STEP for P(drift >= 0 over a window); Obs 334/337 claim a
# refreshed rate kappa' = 0.079/step for D-sojourns. The longest D-capable
# window in a T-step orbit sits at m* ~ log2(#windows)/c: c = 0.2113 -> m* ~ 55,
# c = 0.079 -> m* ~ 146 at T = 3000. Sharply distinguishable.
# D-capable window = consecutive steps with sum(e) <= alpha * length
# (nonnegative log2-drift).
# (i) exact-iid model (e ~ Geom(1/2), proven law of the Haar system), 4000 runs;
# (ii) real integer orbits (1400-bit starts), 300 runs — consistency check.
import numpy as np
from math import log2
import random

ALPHA = log2(3.0)
T = 3000

def longest_window(es):
    # longest j-i with sum_{t in (i, j]} (e_t - alpha) <= 0
    P = np.concatenate([[0.0], np.cumsum(es - ALPHA)])
    runmax = np.maximum.accumulate(P)
    best = 0
    for j in range(1, len(P)):
        # earliest i with runmax[i] >= P[j]
        i = np.searchsorted(runmax, P[j], side='left')
        if i < j and j - i > best:
            best = j - i
    return best

rng = np.random.default_rng(364)
print(f"(i) exact iid model, {T} steps per run:")
ms = []
for _ in range(4000):
    es = rng.geometric(0.5, size=T).astype(float)
    ms.append(longest_window(es))
ms = np.array(ms)
# m* ~ (1/c) * log2(T^2/2) roughly (windows ~ T^2/2 but overlapping: effective T)
# use the standard Erdos-Renyi longest-run law: m*/log2(T) -> 1/c
c_est = log2(T) / ms.mean()
print(f"    longest D-window: mean {ms.mean():.1f}, median {np.median(ms):.0f}, "
      f"p10-p90 {np.percentile(ms,10):.0f}-{np.percentile(ms,90):.0f}")
print(f"    Erdos-Renyi estimate c = log2(T)/E[m*] = {c_est:.4f} bits/step")
print(f"    candidates: I(0) = 0.2113 (Obs 180) | kappa' = 0.079 (Obs 334)")

print()
print(f"(ii) real integer orbits (1400-bit starts), 300 runs:")
ms2 = []
random.seed(364)
for _ in range(300):
    n = random.getrandbits(1400) | 1
    es = np.empty(T)
    for t in range(T):
        n = 3*n + 1
        e = (n & -n).bit_length() - 1
        n >>= e
        es[t] = e
    ms2.append(longest_window(es))
ms2 = np.array(ms2)
c2 = log2(T) / ms2.mean()
print(f"    longest D-window: mean {ms2.mean():.1f}, median {np.median(ms2):.0f}")
print(f"    c = {c2:.4f} bits/step  (integer orbits vs model: "
      f"{'consistent' if abs(c2 - c_est) < 0.03 else 'DEVIATES'})")

# direct tail fit as an independent estimate: P(window >= m) over all runs
print()
allm = np.concatenate([ms, ms2])
qs = np.arange(int(np.percentile(allm, 20)), int(np.percentile(allm, 95)))
frac = [(allm >= m).mean() for m in qs]
sl = np.polyfit(qs, np.log2(np.maximum(frac, 1e-12)), 1)[0]
print(f"direct tail fit of P(m* >= m): slope = {sl:.4f} per step "
      f"(NB: extreme-value tail, slope ~ -c)")
