"""R341-348: THE BACKWARD WAVEFRONT (user's program).
d(n) = total T-steps (n->n/2 | (3n+1)/2) to reach 1. Wavefront B_t = {n: d(n)<=t}.
(1) SHAPE: for each t, min-excluded and fraction reached at scales 2^16..2^24;
    boundary slope c*(x) = d_max(x)/log2(x) -- does the edge stay linear/bounded?
(2) LAGGARDS: record-setting n (d(n) > d(m) for all m<n) -- the numbers that stay
    excluded longest. Their pattern: ratio, binary tail, ternary tail, family (a,k).
(3) Connect to known machinery: are laggards ladder-riders (trailing-ones fuel)?"""
import sys
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")

N = 1 << 24
d = np.full(N, -1, dtype=np.int32)
d[1] = 0
# iterative memo
for n in range(2, N):
    if d[n] >= 0: continue
    path = []
    m = n
    while m >= N or d[m] < 0:
        path.append(m)
        m = m // 2 if m % 2 == 0 else (3 * m + 1) // 2
    base = d[m]
    for x in reversed(path):
        base += 1
        if x < N: d[x] = base

print("(1) wavefront shape: t | min-excluded | frac reached <=2^20 | <=2^24")
tot20 = (1 << 20) - 1; tot24 = N - 1
dmax = int(d.max())
for t in range(10, dmax + 1, 25):
    reached = d[1:] <= t
    # min excluded
    exc = np.flatnonzero(~reached)
    mn = exc[0] + 1 if len(exc) else None
    f20 = reached[:tot20].mean(); f24 = reached.mean()
    print(f"  t={t:>3} | min-excl {mn if mn else '-':>9} | {f20:.4f} | {f24:.4f}")

print(f"\nmax d(n) for n<2^24: {dmax}")

print("\n(2) LAGGARD RECORDS (d(n) sets a new record):")
def tail_ones(x):
    k = 0
    while x & 1: x >>= 1; k += 1
    return k
def tern(x):
    s = ""
    while x: s = str(x % 3) + s; x //= 3
    return s
rec = -1
print(f"{'n':>10} {'d(n)':>5} {'d/log2n':>8} {'bin-tail k':>10} {'family a':>9} {'tern tail':>10}")
laggards = []
for n in range(2, N):
    if d[n] > rec:
        rec = d[n]
        k = tail_ones(n)
        a = (n + 1) >> k if n % 2 else 0
        tt = tern(n)[-6:] if n % 2 else "-"
        laggards.append(n)
        if n > 20:
            print(f"{n:>10} {d[n]:>5} {d[n]/np.log2(n):>8.3f} {k:>10} {a:>9} {tt:>10}")
print("\nlaggard ratios n_{i+1}/n_i:", [f"{b/a:.2f}" for a, b in zip(laggards[-8:], laggards[-7:])])
