"""R433-444: laggard frontier extended to 2^25. New champions, c* growth,
seed-vs-chain classification of the new records."""
import sys
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")

N = 1 << 25
d = np.full(N, -1, dtype=np.int32)
d[1] = 0
for n in range(2, N):
    if d[n] >= 0: continue
    path = []; m = n
    while m >= N or d[m] < 0:
        path.append(m); m = m // 2 if m % 2 == 0 else (3 * m + 1) // 2
    base = d[m]
    for x in reversed(path):
        base += 1
        if x < N: d[x] = base

recs = []
best = -1
for n in range(2, N):
    if d[n] > best: best = d[n]; recs.append(n)
print("records beyond 2^24 (new territory):")
def tail_ones(x):
    k = 0
    while x & 1: x >>= 1; k += 1
    return k
for n in recs:
    if n >= (1 << 24):
        print(f"  n={n:>10} d={d[n]} c*={d[n]/np.log2(n):.3f} k={tail_ones(n)}")
last = recs[-1]
print(f"\nfrontier at 2^25: champion {last}, d={d[last]}, c* = {d[last]/np.log2(last):.3f}")
print("c* sequence at champions:", [f"{d[r]/np.log2(r):.2f}" for r in recs[-6:]])
