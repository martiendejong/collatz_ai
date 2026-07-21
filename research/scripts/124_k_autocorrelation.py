"""
124_k_autocorrelation.py
=========================
K-value autocorrelation along Collatz orbits.

Key question: are consecutive macro-step sizes K (= v2(n+1)) independent,
or do they have memory? Memory would create structure that affects mixing.

Also: compute the Lyapunov exponent (expected log-ratio per step) and
the VARIANCE of the log-ratio to understand orbit fluctuations.
"""
import numpy as np
from collections import defaultdict
import random

def v2(x):
    if x == 0: return 999
    c = 0
    while x % 2 == 0: x //= 2; c += 1
    return c

def macro_step(n):
    K = v2(n+1); m = (n+1) >> K; x = m*(3**K)-1; l = v2(x)
    return x >> l, K, l

def collatz_orbit_K(n_start, max_steps=100000):
    """Run Collatz orbit and collect K values and log-ratios."""
    K_vals = []
    l_vals = []
    log_ratios = []
    n = n_start
    for _ in range(max_steps):
        K = v2(n+1)
        if n == 1: break
        n_out, K_val, l_val = macro_step(n)
        K_vals.append(K_val)
        l_vals.append(l_val)
        log_ratio = K_val * np.log(3) - (K_val + l_val) * np.log(2)
        log_ratios.append(log_ratio)
        n = n_out
    return K_vals, l_vals, log_ratios

print("=" * 70)
print("PART 1: K DISTRIBUTION AND AUTOCORRELATION (multiple starting points)")
print("=" * 70)

# Need very large starting numbers so orbits run MANY macro-steps before reaching 1.
# Expected macro-steps = bits * 1.74. For 5000-bit numbers: ~8700 steps.
import random as _r
_r.seed(42)
starts = [
    _r.getrandbits(5000) | 1,   # 5000-bit random odd
    _r.getrandbits(5000) | 1,
    _r.getrandbits(5000) | 1,
    2**5000 - 1,                 # all ones in binary (5000 bits)
    2**4999 + _r.getrandbits(4999) | 1,
]

all_K = []
all_l = []
all_lr = []
for start in starts:
    Ks, ls, lrs = collatz_orbit_K(start, max_steps=50000)
    all_K.extend(Ks)
    all_l.extend(ls)
    all_lr.extend(lrs)
    print(f"n~2^{start.bit_length()-1}: {len(Ks)} macro-steps, mean K={np.mean(Ks):.4f}, "
          f"mean Lyapunov={np.mean(lrs):.4f}")

K_arr = np.array(all_K, dtype=float)
l_arr = np.array(all_l, dtype=float)
lr_arr = np.array(all_lr, dtype=float)

print(f"\nPooled {len(K_arr)} macro-steps:")
print(f"  E[K] = {K_arr.mean():.5f} (theory: 2.000)")
print(f"  E[l0] = {l_arr.mean():.5f} (theory: 2.000)")
print(f"  Var[K] = {K_arr.var():.5f} (theory: 2.000)")
print(f"  E[Lyapunov] = {lr_arr.mean():.5f}")
print(f"     = E[K]*log3 - (E[K]+E[l0])*log2 = "
      f"{K_arr.mean()*np.log(3) - (K_arr.mean()+l_arr.mean())*np.log(2):.5f} (predicted)")
print(f"  Std[Lyapunov] = {lr_arr.std():.5f}")

print("\nK distribution:")
K_counts = defaultdict(int)
for k in K_arr.astype(int): K_counts[k] += 1
total = len(K_arr)
for k in sorted(K_counts)[:12]:
    freq = K_counts[k]/total
    theory = 1/2**k if k >= 1 else 0
    print(f"  K={k:2d}: freq={freq:.5f} theory=1/2^{k}={theory:.5f} ratio={freq/theory:.3f}")

print("\nK autocorrelation:")
K_mean = K_arr.mean(); K_var = K_arr.var()
print(f"  {'Lag':>5} {'ACF':>12} {'ACF_l0':>12}")
for lag in range(1, 15):
    if lag >= len(K_arr): break
    acf_K = np.corrcoef(K_arr[:-lag], K_arr[lag:])[0,1]
    acf_l = np.corrcoef(l_arr[:-lag], l_arr[lag:])[0,1]
    print(f"  {lag:>5} {acf_K:>12.6f} {acf_l:>12.6f}")

print("\nCross-correlation K_t vs l0_{t+1} (do high-K steps predict l0?):")
for lag in range(0, 5):
    if lag >= len(K_arr): break
    if lag == 0:
        cc = np.corrcoef(K_arr, l_arr)[0,1]
    else:
        cc = np.corrcoef(K_arr[:-lag], l_arr[lag:])[0,1]
    print(f"  lag {lag}: corr(K_t, l0_{{t+{lag}}}) = {cc:.6f}")

print()
print("=" * 70)
print("PART 2: JOINT DISTRIBUTION OF (K, l0)")
print("=" * 70)
print("\nJoint distribution P(K=k, l0=l) vs independence P(K=k)*P(l0=l):")
joint = defaultdict(int)
for k, l in zip(K_arr.astype(int), l_arr.astype(int)):
    joint[(k,l)] += 1
total = len(K_arr)

# Show entries with freq > 0.5%
print(f"  {'K':>4} {'l0':>4} {'Observed':>12} {'Independent':>12} {'Ratio':>8}")
for (k, l), count in sorted(joint.items(), key=lambda x: -x[1])[:20]:
    obs = count/total
    ind = (K_counts[k]/total) * (sum(1 for li in l_arr.astype(int) if li==l)/total)
    if obs > 0.001:
        print(f"  {k:>4} {l:>4} {obs:>12.5f} {ind:>12.5f} {obs/ind:>8.3f}")

print()
print("=" * 70)
print("PART 3: LYAPUNOV EXPONENT ANALYSIS")
print("=" * 70)
print(f"\nLyapunov exponent = E[K]*log(3) - (E[K]+E[l0])*log(2)")
print(f"  If = 0: orbits neutral (boundary case)")
print(f"  If < 0: orbits contract (supports Collatz)")
print(f"  Measured: {lr_arr.mean():.6f}")
print(f"  Theoretical: 2*log(3) - 4*log(2) = {2*np.log(3) - 4*np.log(2):.6f}")
print(f"  Typical fluctuation per step: std = {lr_arr.std():.4f}")
print(f"  Steps to shrink to 1 from n: ~log(n)/|Lyapunov| = {np.log(2**60)/abs(lr_arr.mean()):.1f} for n~2^60")
print()

# Distribution of log-ratio
print("Log-ratio distribution (K*log3 - (K+l0)*log2):")
bins = np.linspace(-3, 3, 25)
hist, _ = np.histogram(lr_arr, bins=bins)
for i in range(len(hist)):
    bar = '#' * (hist[i] * 40 // max(hist))
    print(f"  [{bins[i]:+.2f},{bins[i+1]:+.2f}): {bar} ({100*hist[i]/total:.1f}%)")
