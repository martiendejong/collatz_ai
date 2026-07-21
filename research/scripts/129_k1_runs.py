"""
129_k1_runs.py
===============
Structure of consecutive K=1 runs in Collatz orbits.

After a K>=2 macro-step, how long is the subsequent run of K=1 steps?
Does the run length depend on K (the preceding step size)?

Also: what is the JOINT distribution of (K_prev, run_length)?

Context: From orbit of n=127: after K=7 step -> 1093, then K=1 -> 205,
then many K=1 steps: 205->77->29->11 (all K=1). This suggests K-1 runs
may have structure related to 3-adic structure.

Key: K=1 steps are the "boring" Collatz steps (n->3n/2 when n=odd and (n+1)/2 odd).
     K=2+ steps are the "rich" steps (multiple consecutive halvings).
"""
import numpy as np
from collections import defaultdict
import random as _r

def v2(x):
    if x == 0: return 999
    c = 0
    while x % 2 == 0: x //= 2; c += 1
    return c

def macro_step(n):
    K = v2(n+1); m = (n+1) >> K; x = m*(3**K)-1; l = v2(x)
    return x >> l, K, l

print("=" * 70)
print("PART 1: K=1 RUN LENGTH DISTRIBUTION")
print("=" * 70)
print()

_r.seed(42)
n = _r.getrandbits(3000) | 1  # 3000-bit random odd number

# Collect K values and identify K=1 runs
K_sequence = []
for _ in range(100000):
    n_out, K_val, l_val = macro_step(n)
    K_sequence.append((K_val, l_val))
    n = n_out
    if n < 2:
        n = _r.getrandbits(3000) | 1

# Find all maximal K=1 runs
run_lengths = []
current_run = 0
for K, l in K_sequence:
    if K == 1:
        current_run += 1
    else:
        if current_run > 0:
            run_lengths.append(current_run)
        current_run = 0
if current_run > 0:
    run_lengths.append(current_run)

run_dist = defaultdict(int)
for r in run_lengths:
    run_dist[r] += 1

total_runs = len(run_lengths)
total_K1_steps = sum(run_lengths)
total_steps = len(K_sequence)

print(f"Total macro-steps: {total_steps}")
print(f"Total K=1 steps: {total_K1_steps} ({100*total_K1_steps/total_steps:.1f}%)")
print(f"Total K>=2 steps: {total_steps-total_K1_steps} ({100*(total_steps-total_K1_steps)/total_steps:.1f}%)")
print(f"Total K=1 runs: {total_runs}")
print(f"Mean K=1 run length: {np.mean(run_lengths):.4f}")
print()

print("K=1 run length distribution:")
print(f"  {'Run len':>8} {'Count':>8} {'Prob':>10} {'CDF':>8} {'Geom(p)':>10}")
cdf = 0
# Estimate geometric parameter
K1_frac = 0.5  # P(K=1) = 1/2
# P(run >= r) = P(r consecutive K=1 steps) = (1/2)^r
# P(run = r) = P(K=1)^r * P(K>=2) = (1/2)^r * (1/2) -- for geometric
# But runs are between K>=2 steps, so P(run=r) = P(K=1)^r * P(K>=2) = (1/2)^r * (1/2)
for rl in sorted(run_dist)[:25]:
    prob = run_dist[rl]/total_runs
    cdf += prob
    geom_pred = (0.5**rl) * 0.5  # geometric with p=1/2
    print(f"  {rl:>8} {run_dist[rl]:>8} {prob:>10.5f} {cdf:>8.4f} {geom_pred:>10.5f}")

print()
print("=" * 70)
print("PART 2: WHAT K>=2 VALUE PRECEDES EACH K=1 RUN?")
print("=" * 70)
print()

# For each K=1 run, what was the K>=2 step that preceded it?
# Also: what K>=2 step follows the run?
run_contexts = []  # (K_before, run_length, K_after)
prev_K = None
prev_was_high = False
temp_run = 0
first_high_K = None

for K, l in K_sequence:
    if K == 1:
        if prev_was_high:
            temp_run = 1
            first_high_K = prev_K
        elif temp_run > 0:
            temp_run += 1
        else:
            temp_run = 1
            first_high_K = None  # run at start
        prev_was_high = False
    else:
        if temp_run > 0 and first_high_K is not None:
            run_contexts.append((first_high_K, temp_run, K))
        temp_run = 0
        first_high_K = None
        prev_was_high = True
    prev_K = K

# Distribution of K_before for different run lengths
print("Mean K=1 run length conditioned on K_before:")
Kbefore_groups = defaultdict(list)
for K_before, run_len, K_after in run_contexts:
    Kbefore_groups[K_before].append(run_len)

print(f"  {'K_before':>10} {'Count':>8} {'MeanRun':>10} {'Var':>8} {'P(run=1)':>10}")
for K in sorted(Kbefore_groups)[:12]:
    runs = Kbefore_groups[K]
    print(f"  {K:>10} {len(runs):>8} {np.mean(runs):>10.4f} {np.var(runs):>8.4f} {(np.array(runs)==1).mean():>10.4f}")

print()
print("K=1 run length distribution conditioned on K_before=1,2,3,4,5+:")
for K_thr in [1, 2, 3, 4, 5]:
    if K_thr == 5:
        runs = [r for (Kb,r,_) in run_contexts if Kb >= 5]
        label = "K>=5"
    else:
        runs = [r for (Kb,r,_) in run_contexts if Kb == K_thr]
        label = f"K={K_thr}"
    if not runs: continue
    run_arr = np.array(runs)
    print(f"  {label}: n={len(runs)}, mean={run_arr.mean():.3f}, "
          f"P(1)={( run_arr==1).mean():.3f}, P(2)={( run_arr==2).mean():.3f}, "
          f"P(3)={( run_arr==3).mean():.3f}, P(>=4)={(run_arr>=4).mean():.3f}")

print()
print("=" * 70)
print("PART 3: K=1 AUTOCORRELATION WITHIN RUNS")
print("=" * 70)
print()
print("Within a K=1 run, are consecutive l0 values correlated?")
# Collect l0 values within K=1 runs
l0_in_runs = []
current_l0_run = []
for K, l in K_sequence:
    if K == 1:
        current_l0_run.append(l)
    else:
        if len(current_l0_run) >= 2:
            l0_in_runs.append(current_l0_run[:])
        current_l0_run = []

# Pool all within-run l0 pairs (consecutive l0 values)
l0_t = []; l0_t1 = []
for run in l0_in_runs:
    for i in range(len(run)-1):
        l0_t.append(run[i]); l0_t1.append(run[i+1])

if l0_t:
    l0_t = np.array(l0_t); l0_t1 = np.array(l0_t1)
    corr = np.corrcoef(l0_t, l0_t1)[0,1]
    print(f"Consecutive l0 autocorrelation within K=1 runs: {corr:.6f}")

    # Also: l0 distribution within K=1 runs (should still be geometric)
    l0_all_in_runs = []
    for run in l0_in_runs:
        l0_all_in_runs.extend(run)
    l0_arr = np.array(l0_all_in_runs)
    print(f"\nl0 distribution within K=1 runs:")
    from collections import Counter
    cnt = Counter(l0_arr)
    total_l0 = len(l0_arr)
    for l in sorted(cnt)[:8]:
        obs = cnt[l]/total_l0
        theory = 1/2**l
        print(f"  l0={l}: {obs:.5f} (theory 1/2^l={theory:.5f}, ratio={obs/theory:.3f})")

print()
print("=" * 70)
print("PART 4: EFFECTIVE STEP SIZE IN K=1 RUNS")
print("=" * 70)
print()
print("A K=1, l0=2 step: n -> (3(n+1)/2 - 1)/4 = (3n/2 + 3/2 - 1)/4")
print("  Actually: n -> 3m - 1 where m=(n+1)/2, then /2^l0")
print("  For K=1, l0=1: n_out = (3(n+1)/2 - 1)/2 = (3n+1)/4")
print("  For K=1, l0=2: n_out = (3(n+1)/2 - 1)/4 = (3n+1)/8")
print("  For K=1, l0=r: n_out = (3n+1) / 2^{r+1}")
print()
print("The l0 distribution for K=1 is: P(l0=r) = 1/2^r (geometric).")
print("Mean l0 for K=1: E[l0] = 2.")
print("Effective multiplier for K=1 step: 3 / 2^{l0+1}")
print("  l0=1: 3/4 = 0.75 (CONTRACT)")
print("  l0=2: 3/8 = 0.375 (CONTRACT fast)")
print("  l0=3: 3/16 = 0.1875 (CONTRACT very fast)")
print()
print("Lyapunov for K=1: E[log ratio] = log(3) - (1+E[l0])*log(2) = log(3)-3*log(2)")
print(f"  = {__import__('math').log(3) - 3*__import__('math').log(2):.6f}")
print("  (negative: K=1 steps always contract)")
print()

# Empirical: after a run of r K=1 steps, what is the net log ratio?
net_ratios = []
for run_contexts_item in [(Kb,r,Ka) for (Kb,r,Ka) in run_contexts if r >= 1]:
    Kb, rl, Ka = run_contexts_item
    # Find the l0 values in this run
    pass  # Would need to store l0 values separately

print("Net Lyapunov per K=1 step (from K sequence):")
K1_l0_vals = [l for K,l in K_sequence if K==1]
K1_lyap = [__import__('math').log(3) - (1+l)*__import__('math').log(2) for l in K1_l0_vals]
print(f"  Mean: {__import__('numpy').mean(K1_lyap):.6f}")
print(f"  Std:  {__import__('numpy').std(K1_lyap):.6f}")
print(f"  Theory: log(3)-3*log(2) = {__import__('math').log(3)-3*__import__('math').log(2):.6f}")
