"""
235_cycle_diophantine.py
========================
SP1A: Near-miss (k,h) pairs where 2^h is closest to 3^k from above.
SP1B: Lower bound on n0 — which k-values are ruled out by computational bound 2^68.

Cycle equation for a Collatz cycle with k odd steps and h total halvings:
  n0 * (2^h - 3^k) = S(h_1,...,h_k)

where S = sum_{j=0}^{k-1} 3^j * 2^{H_j}, H_j = h_{j+1}+...+h_k (back-loaded halvings).
Valid cycles require: 2^h > 3^k, n0 odd positive integer, h_i >= 1, sum h_i = h.

SP1A output: for each k=1..200, the unique h_k = ceil(k*log2(3)) with eps_k = h_k - k*log2(3).
Near-miss pairs (small eps) come from convergents of CF of log2(3).

SP1B output: lower bound on n0 from the cycle equation combined with S_min.
The minimum S over all halving patterns (h_i>=1, sum=h) is achieved by
front-loading: h_1 = h-k+1, h_2=...=h_k=1. Then:
  S_min = 2^h + sum_{j=1}^{k-1} 3^j * 2^{k-j} = 2^h + 2*3^k - 3*2^k.
"""
import sys
from math import log, log2, ceil, floor

try:
    import mpmath
    mpmath.mp.dps = 60
    LOG2_3 = mpmath.log(3, 2)
    HAS_MPMATH = True
except ImportError:
    LOG2_3 = log2(3)  # fallback; less accurate for large k
    HAS_MPMATH = False

print("235: Cycle Diophantine constraints (SP1A + SP1B)")
print(f"log2(3) = {float(LOG2_3):.15f}")
if HAS_MPMATH:
    print("  (mpmath 60-digit precision)")
print("=" * 72)
sys.stdout.flush()

# --- CF convergents of log2(3) ---
# Compute via the standard recursive algorithm on the fractional part of log2(3)
# Known CF: [1; 1, 1, 2, 2, 3, 1, 5, 2, 23, 2, 2, 1, 1, 55, 1, 4, 3, 1, 1, ...]
# Convergents h_n/k_n with h_n/k_n alternately above/below log2(3).

def compute_cf_convergents(target, max_k=2000):
    """
    Convergents (h_n, k_n) of target where h_n/k_n -> target.
    Returns list of (k_n, h_n, sign) where sign=+1 means h_n/k_n > target.
    """
    x = float(target)
    p_prev, p_curr = 1, int(x)
    q_prev, q_curr = 0, 1
    convergents = [(q_curr, p_curr)]
    x_rem = x - int(x)
    for _ in range(200):
        if x_rem < 1e-14:
            break
        x_rem = 1.0 / x_rem
        a = int(x_rem)
        x_rem -= a
        p_prev, p_curr = p_curr, a * p_curr + p_prev
        q_prev, q_curr = q_curr, a * q_curr + q_prev
        if q_curr > max_k:
            break
        convergents.append((q_curr, p_curr))
    return convergents

convergents = compute_cf_convergents(float(LOG2_3), max_k=2000)
conv_set = set((k, h) for k, h in convergents)

print("\n--- Continued fraction convergents of log2(3) ---")
print(f"  {'n':>3}  {'k':>6}  {'h':>7}  {'h/k':>12}  {'eps=h-k*log2(3)':>18}  {'above?'}")
for i, (k, h) in enumerate(convergents):
    if HAS_MPMATH:
        lk = float(mpmath.mpf(k) * LOG2_3)
    else:
        lk = k * float(LOG2_3)
    eps = h - lk
    above = "YES" if eps > 0 else "no"
    print(f"  {i:>3}  {k:>6}  {h:>7}  {h/k:>12.9f}  {eps:>18.9f}  {above}")
    if i >= 18:
        print("  ...")
        break
print()
sys.stdout.flush()

# --- SP1A: For each k=1..200, the unique h_k with 2^h > 3^k smallest ---
print("--- SP1A: Near-miss (k,h) pairs with smallest eps = h - k*log2(3) ---")
print(f"  {'k':>5}  {'h':>6}  {'eps':>12}  {'log2(gap)':>12}  {'conv?':>6}")

results = []
for k in range(1, 201):
    if HAS_MPMATH:
        lk = mpmath.mpf(k) * LOG2_3
        h = int(mpmath.ceil(lk))
        if h <= float(lk) - 0.5:  # guard: ceil can be equal for exact values
            h += 1
        eps = float(h - lk)
        # log2(gap) = log2(2^h - 3^k) = log2(3^k * (2^eps - 1))
        #           = k*log2(3) + log2(2^eps - 1)
        if eps < 1e-10:
            log2_gap = float(k * LOG2_3) + log2(eps * log(2) + 1e-300)
        else:
            log2_gap = float(k * LOG2_3) + log2(2**eps - 1)
    else:
        lk = k * float(LOG2_3)
        h = ceil(lk)
        if h == lk:
            h += 1
        eps = h - lk
        log2_gap = lk + log2(max(2**eps - 1, 1e-300))

    is_conv = (k, h) in conv_set
    results.append((k, h, eps, log2_gap, is_conv))

# Print top-20 smallest eps (most dangerous near-misses)
sorted_eps = sorted(results, key=lambda x: x[2])
print("\n  Top-20 pairs with smallest eps (most constrained n0):")
for k, h, eps, lg, is_conv in sorted_eps[:20]:
    tag = "*CONV*" if is_conv else ""
    print(f"  k={k:>5}  h={h:>6}  eps={eps:>12.8f}  log2(gap)={lg:>9.3f}  {tag}")

print()
sys.stdout.flush()

# --- SP1B: Lower bound on n0 from S_min ---
# S_min formula (front-loaded halvings):
#   S_min = 2^h + 2*3^k - 3*2^k
# For large k: S_min ≈ 2^h + 2*3^k ≈ 3*3^k (since 2^h ≈ 3^k)
# n0_min = S_min / (2^h - 3^k) = S_min / gap
# In log2 units:
#   log2(n0_min) = log2(S_min) - log2(gap)
#                ≈ log2(3*3^k) - log2(gap)
#                = log2(3) + k*log2(3) - log2(gap)

print("--- SP1B: Lower bound on n0 = S_min / gap ---")
print("  S_min (front-loaded halvings) = 2^h + 2*3^k - 3*2^k")
print("  log2(n0_min) = log2(S_min) - log2(gap)")
print()
print(f"  Computational exclusion bound: n0 <= 2^68 = {2**68:.3e}")
print(f"  If log2(n0_min) > 68: this (k,h) is outside the computational range.")
print()

TARGET_LOG2 = 68.0

print(f"  {'k':>5}  {'h':>6}  {'eps':>10}  {'log2(S_min)':>12}  {'log2(gap)':>11}  {'log2(n0_min)':>13}  {'>68?'}")

K_CRIT = None
for k, h, eps, log2_gap, is_conv in results[:200]:
    # log2(S_min) = log2(2^h + 2*3^k - 3*2^k)
    # For k >= 2: 2^h ≈ 3^k, 3*2^k << 2*3^k
    # log2(S_min) ≈ log2(3) + k*log2(3)  [dominant term is 2*3^k ≈ 2*2^h/1 ≈ 2^h]
    # More carefully: S_min/3^k = 2^{h-k*log2(3)} + 2 - 3*(2/3)^k
    #                           = 2^eps + 2 - 3*(2/3)^k
    #                           ≈ 1 + 2 = 3  for large k (eps small, (2/3)^k tiny)
    kl3 = k * float(LOG2_3)
    # log2(S_min) = log2(2^h + 2*3^k - 3*2^k)
    # Use: log2(3^k * (2^eps + 2 - 3*(2/3)^k))
    ratio = 2**eps + 2.0 - 3.0 * (2.0/3.0)**k
    if ratio <= 0:
        ratio = 1e-10  # shouldn't happen for valid k
    log2_Smin = kl3 + log2(ratio)

    log2_n0_min = log2_Smin - log2_gap

    above = "YES" if log2_n0_min > TARGET_LOG2 else "no"
    if log2_n0_min > TARGET_LOG2 and K_CRIT is None:
        K_CRIT = k

    if k <= 30 or is_conv or k % 20 == 0:
        tag = "*" if is_conv else ""
        print(f"  k={k:>4}  h={h:>5}  eps={eps:>10.6f}  {log2_Smin:>12.2f}  {log2_gap:>11.2f}  {log2_n0_min:>13.2f}  {above} {tag}")

print()
print(f"  First k where n0_min > 2^68: k = {K_CRIT}")
print()
print("Note: S_min bound is tight only for front-loaded halvings (h_1=h-k+1, rest=1).")
print("For generic halving patterns, n0 may be larger. This is a *lower* bound on n0_min.")
print()
print("Known result (Simons-de Weger 2003): no non-trivial cycles with k < 35000.")
print("Combined with computational verification (n0 <= 2^68): no non-trivial cycles")
print("exist for k < K_CRIT OR n0 <= 2^68.")

# --- Summary table for paper ---
print()
print("--- Summary table (for paper): convergent pairs k<=200 ---")
print(f"  {'k':>5}  {'h':>6}  {'eps':>10}  {'log2(n0_min)':>13}  note")
for k, h, eps, log2_gap, is_conv in results:
    if is_conv:
        kl3 = k * float(LOG2_3)
        ratio = 2**eps + 2.0 - 3.0 * (2.0/3.0)**k
        if ratio <= 0: ratio = 1e-10
        log2_Smin = kl3 + log2(ratio)
        log2_n0_min = log2_Smin - log2_gap
        note = "convergent"
        print(f"  k={k:>4}  h={h:>5}  eps={eps:>10.6f}  {log2_n0_min:>13.2f}  {note}")

print()
print("done")
