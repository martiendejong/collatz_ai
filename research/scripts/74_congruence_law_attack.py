"""R1516-1540: ATTACKING THE CRITICAL-WINDOW CONGRUENCE LAW.
(1) CLEARANCE: min over all words of distance |W mod D| to 0 (cyclic), per r.
    The quantitative form of the law: how close does any word come to a
    multiple of D? Meet-in-the-middle + sorted nearest search, r <= 24.
(2) NON-CRITICAL windows j = ceil+1: is 0 reachable there? (is criticality
    essential, or is the law about all windows?)
(3) r=18 FULL-D DP via bigint bitsets (D = 149,450,423).
"""
import sys, math, bisect
from math import comb
from collections import defaultdict
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
ALPHA = math.log2(3)

# ---------- (1) clearance ----------
def clearance(r, j):
    D = 2 ** j - 3 ** r
    r1 = r // 2; r2 = r - r1
    def half(rr, smax):
        out = defaultdict(list)
        def rec(parts_left, s, J, val):
            if parts_left == 0:
                out[s].append(val % D); return
            i = rr - parts_left
            term = pow(3, rr - 1 - i, D) * pow(2, J, D)
            for js in range(1, smax - s - (parts_left - 1) + 1):
                rec(parts_left - 1, s + js, J + js, val + term)
        rec(rr, 0, 0, 0)
        return out
    L = half(r1, j - r2); R = half(r2, j - r1)
    p3 = pow(3, r2, D)
    best = D
    for s, la in L.items():
        need = j - s
        rb = R.get(need)
        if not rb: continue
        p2s = pow(2, s, D)
        tb = np.sort((np.array(rb, dtype=object) * p2s) % D).astype(object)
        lav = (np.array(la, dtype=object) * p3) % D
        for a in lav:
            t = (D - a) % D          # want tb close to -a mod D
            i = bisect.bisect_left(tb, t)
            for cand in (tb[i % len(tb)], tb[(i - 1) % len(tb)]):
                x = (a + cand) % D
                best = min(best, x, D - x)
    return D, int(best)

print("(1) CLEARANCE per critical window: min cyclic distance of any word value to 0 mod D")
print(f"{'r':>3s} {'j':>3s} {'D':>15s} {'clearance':>12s} {'clr/D':>9s} {'log2 clr':>9s}")
for r in range(3, 21):
    j = math.ceil(r * ALPHA)
    if 2 ** j - 3 ** r <= 0: j += 1
    D, c = clearance(r, j)
    print(f"{r:3d} {j:3d} {D:15d} {c:12d} {c/D:9.5f} {math.log2(c) if c>0 else -1:9.2f}")

# ---------- (2) non-critical windows ----------
def reach0(r, j, D):
    dp = np.zeros((j, D), dtype=bool)
    dp[0][pow(3, r - 1, D)] = True
    for i in range(1, r):
        coef = pow(3, r - 1 - i, D)
        dp2 = np.zeros((j, D), dtype=bool)
        pref = np.zeros(D, dtype=bool)
        for Jn in range(1, j):
            pref |= dp[Jn - 1]
            if pref.any():
                dp2[Jn] = np.roll(pref, (coef * pow(2, Jn, D)) % D)
        dp = dp2
    out = dp.any(axis=0)
    return bool(out[0])

print("\n(2) NON-CRITICAL windows j = ceil+1: 0 reachable mod D?")
for r in range(3, 13):
    j = math.ceil(r * ALPHA) + 1
    D = 2 ** j - 3 ** r
    ok = reach0(r, j, D)
    print(f"  r={r:2d} j={j:2d} D={D:8d}: 0 {'REACHABLE' if ok else 'BLOCKED'}")

# ---------- (3) r=18 full-D bigint DP ----------
print("\n(3) r=18 critical window, D = 149,450,423 (bigint bitset DP):", flush=True)
r = 18; j = 29; D = 2 ** j - 3 ** r
mask = (1 << D) - 1
def rot(x, c):
    return ((x << c) | (x >> (D - c))) & mask
dp = [0] * j
dp[0] = 1 << (pow(3, r - 1, D) % D)
for i in range(1, r):
    coef = pow(3, r - 1 - i, D)
    dp2 = [0] * j
    pref = 0
    for Jn in range(1, j):
        pref |= dp[Jn - 1]
        if pref:
            dp2[Jn] = rot(pref, (coef * pow(2, Jn, D)) % D)
    dp = dp2
    print(f"  step {i}/{r-1} done", flush=True)
out = 0
for J in range(j): out |= dp[J]
print(f"r=18: 0 mod D {'REACHABLE' if out & 1 else 'BLOCKED'}; reach fraction = {bin(out).count('1')/D:.4f}")
