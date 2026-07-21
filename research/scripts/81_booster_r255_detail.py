"""
81_booster_r255_detail.py
==========================
Detailed transition analysis for r=255 (k=8), the dominant booster.
Focus: distribution of hop lengths, high-k hop probability, exact short-return analysis.
"""
import math, sys, time
from collections import Counter, defaultdict

sys.stdout.reconfigure(encoding='utf-8')

def v2(x):
    c = 0
    while x % 2 == 0: x //= 2; c += 1
    return c

def macro_step(n):
    k = v2(n+1)
    m = (n+1) >> k
    x = m * (3**k) - 1
    l = v2(x)
    return x >> l, k, l

BSet = {27, 55, 63, 83, 95, 103, 127, 159, 169, 191, 207, 223, 239, 253, 255}
DK_THRESHOLD = 3.419

# 1. Exact output distribution from r=255 for first macro-step
# n = 256m - 1 where m is odd (for k=8 exactly: n+1 = 256m, m odd)
# Actually n equiv 255 mod 512 gives exactly k=8
print("=== EXACT 1-STEP OUTPUT DISTRIBUTION FROM r=255 (k=8) ===")
print("n = 256m - 1, m odd, stepping m = 1, 3, 5, ..., 511 (one full period in m mod 256)")
print()

counts_1step = Counter()
period_m = 512  # 256 odd values (m = 1, 3, ..., 511)
BSet_hits_1step = []

for m in range(1, period_m, 2):   # m = 1, 3, ..., 511
    n = 256 * m - 1                # n equiv 255 mod 512, k=8 exactly
    out, k, l = macro_step(n)
    assert k == 8, f"k mismatch at m={m}"
    r_out = out % 256
    counts_1step[r_out] += 1
    if r_out in BSet:
        BSet_hits_1step.append((m, r_out, out, l, k))

total_1step = sum(counts_1step.values())
bset_1step = sum(counts_1step[r] for r in BSet)
print(f"Period: {total_1step} values of m (odd, mod 512)")
print(f"Direct BSet hits (h=1): {bset_1step}/{total_1step} = {100*bset_1step/total_1step:.2f}%")
print()
print("BSet targets from 1-step:")
bset_detail = {r: counts_1step[r] for r in BSet if counts_1step[r] > 0}
for r2, cnt in sorted(bset_detail.items()):
    print(f"  r'={r2:3d}: count={cnt}/{total_1step} = {100*cnt/total_1step:.2f}%  k={v2(r2+1)}")
print()
print("All 1-step BSet transitions (m values):")
for m, r2, out, l, k in sorted(BSet_hits_1step, key=lambda x: x[1]):
    print(f"  m={m:3d}  r'={r2:3d}  output={out}  l={l}")

# 2. The "optimal 2-cycle" r=255 -> r' -> r=255 analysis
print()
print("=== OPTIMAL 2-CYCLES: r=255 -> r' -> r=255 ===")
print("For each r' reachable from r=255 in 1 step,")
print("find m-values where the second hop ALSO lands in BSet")

# For each BSet_hit, compute the 2nd hop
print()
print(f"{'r':>5} {'r_prime':>8} {'r_prime2':>10} {'h_total':>8} {'k_total':>9} {'avg_k':>7}")
print("-"*55)

two_hop_results = []
for m1, r_prime, out1, l1, k1 in BSet_hits_1step:
    n2 = out1
    out2, k2, l2 = macro_step(n2)
    r_pp = out2 % 256
    h_total = 2
    k_total = k1 + k2
    avg_k = k_total / h_total
    two_hop_results.append((r_prime, r_pp, h_total, k_total, avg_k, m1))
    print(f"  r=255 -> r'={r_prime:3d} -> r''={r_pp:3d}  h={h_total}  k_sum={k_total}  avg_k={avg_k:.3f}")

print()

# 3. Extended booster simulation with detailed hop histogram for r=255
print("=== HOP LENGTH HISTOGRAM FOR r=255 (10000 samples) ===")
BASE = 512 * 1_000_000
step = 512  # for k=8 exactly

n_base = (BASE // step) * step + 255
if n_base < BASE:
    n_base += step
_, k_check, _ = macro_step(n_base)
assert k_check == 8

N = 10000
hop_lengths = Counter()
next_booster = Counter()
hop_by_dest = defaultdict(list)  # r' -> list of (steps, k_sum)

converged = 0
for j in range(N):
    n = n_base + j * step
    hops = 0
    ks = 0
    n_cur = n
    while True:
        n_cur, k, l = macro_step(n_cur)
        hops += 1
        ks += k
        if n_cur <= 1:
            converged += 1
            break
        r2 = n_cur % 256
        if r2 in BSet:
            hop_lengths[hops] += 1
            next_booster[r2] += 1
            hop_by_dest[r2].append((hops, ks))
            break
        if hops > 5000:
            break

total_counted = sum(hop_lengths.values())
print(f"N={N}, reached booster={total_counted}, converged={converged}")
print()
print("Hop length distribution:")
for h in sorted(hop_lengths.keys())[:20]:
    cnt = hop_lengths[h]
    bar = '#' * int(50 * cnt / max(hop_lengths.values()))
    print(f"  h={h:3d}: {cnt:5d} ({100*cnt/total_counted:.1f}%) {bar}")
print(f"  ...")

print()
print("Top destination boosters from r=255:")
for r2, cnt in sorted(next_booster.items(), key=lambda x: -x[1]):
    if r2 in hop_by_dest and len(hop_by_dest[r2]) > 0:
        data = hop_by_dest[r2]
        avg_h = sum(h for h,_ in data) / len(data)
        avg_k = sum(ks for _,ks in data) / sum(h for h,_ in data)
        prob = cnt / total_counted
        print(f"  r'={r2:3d} (k={v2(r2+1)}): p={prob:.3f}  avg_h={avg_h:.2f}  avg_k/step={avg_k:.4f}  n={cnt}")

# 4. Best 2-step and 3-step paths from r=255 (h=1+h' or h=2+h')
print()
print("=== BEST COMPOSITE PATHS: r=255 -> r' -> r=255 ===")
print("(using empirical avg transitions)")

# For each r' reachable directly (h=1) from r=255, compute composite avg_k
# if r' then transitions back to r=255
bset_direct_from_255 = {r2: (hop_by_dest[r2] if hop_by_dest[r2] else []) for r2 in BSet if r2 in next_booster}

for r2 in sorted(bset_direct_from_255.keys()):
    # Count h=1 hops from r=255 to r2
    h1_count = sum(1 for h, ks in hop_by_dest.get(r2, []) if h == 1)
    total_r2 = len(hop_by_dest.get(r2, []))
    if h1_count == 0:
        continue
    # For r' -> r=255: need hop data from r2 back to r=255
    # (we don't have this directly, but can estimate from the transition data)
    print(f"  r=255 -h=1-> r={r2:3d} (k={v2(r2+1)}): {h1_count}/{total_r2} hops are h=1")
    # k_sum for 1-step hop: k=8 (from 255) + k_r2 (from r2)
    k_r2 = v2(r2+1)
    print(f"    Combined k for h=2 path: {8+k_r2} -> avg_k={( 8+k_r2)/2:.3f}")

print()
print("=== SUSTAINED HIGH-K CHAIN ANALYSIS ===")
# Find the MAXIMUM avg_k achievable in any finite window of hops
# that starts at r=255 and ends at a booster
print("Looking for windows of consecutive hops with avg_k >= threshold...")

# Collect all hop sequences starting from r=255 (multi-hop)
windows_above = 0
max_window_k = 0.0
max_window_n = 0

for j in range(2000):
    n = n_base + j * step
    # Collect up to 10 hops
    hops_seq = []
    n_cur = n
    for _ in range(10):
        n_cur, k, l = macro_step(n_cur)
        hops_seq.append(k)
        if n_cur <= 1:
            break
        r2 = n_cur % 256
        if r2 in BSet:
            hops_seq = hops_seq  # stop here
            break

    # Check all windows of hops
    for start in range(len(hops_seq)):
        for end in range(start+1, len(hops_seq)+1):
            window = hops_seq[start:end]
            ak = sum(window) / len(window)
            if ak > max_window_k:
                max_window_k = ak
                max_window_n = n
            if ak >= DK_THRESHOLD:
                windows_above += 1

print(f"Total windows checked from {2000} starts: max_avg_k in any window={max_window_k:.4f}")
print(f"Windows with avg_k >= {DK_THRESHOLD}: {windows_above}")

print()
print(f"[Done]")
