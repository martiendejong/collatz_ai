"""R266-271: SPINE ANATOMY & CORRIDOR ENTROPY.
(1) Is the depth-d spine node simply the minimum odd value at backward depth d?
    Build exact backward layers (min-value per depth) and compare with the
    empirical top-1 of R255.
(2) Layer entropy H_d of the visit distribution: growth rate = corridor entropy
    (bits/layer); compare with log2 of branching (2 preimages/window minus springs).
(3) Spine growth factor: value ratio of consecutive spine nodes."""
import sys, random, math
from collections import Counter
sys.stdout.reconfigure(encoding="utf-8")
random.seed(266)

# --- exact backward layers: depth = number of odd predecessors steps from 1
# preimages of odd v (excluding v itself): n = (v*2^w - 1)/3 odd, w >= 1
def preimages(v, cap=10**12):
    out = []
    w = 1
    while v * (1 << w) <= 3*cap:
        x = v * (1 << w) - 1
        if x % 3 == 0:
            n = x // 3
            if n % 2 == 1 and n > 1: out.append(n)
        w += 1
    return out

# BFS by depth, tracking min value per depth (beam on smallest values suffices:
# min at depth d+1 is a preimage of SOME node at depth d; keep the 2000 smallest per layer)
layer = {5}
min_per_depth = []
for d in range(1, 15):
    nxt = set()
    for v in layer:
        for p in preimages(v, cap=10**10):
            if p % 3 == 0: continue  # springs have no further ancestry but ARE nodes
            nxt.add(p)
        # include springs as nodes (they occupy the layer even if childless)
        for p in preimages(v, cap=10**10):
            nxt.add(p)
    layer = set(sorted(nxt)[:4000])
    min_per_depth.append(min(nxt) if nxt else None)

print("exact minimum odd value per backward depth (from 5 = depth 1):")
print("  depth:", list(range(2, 16)))
print("  min:  ", min_per_depth)

# --- empirical spine (top-visited per depth) for comparison
N = 120_000
layers = [Counter() for _ in range(15)]
for _ in range(N):
    n = random.randrange(3, 1 << 44) | 1
    odds = [n]; m = n
    while m != 1:
        m = 3*m+1 if m % 2 else m//2
        if m % 2 and m != 1: odds.append(m)
    for d in range(1, 16):
        if len(odds) >= d: layers[d-1][odds[-d]] += 1

print("\nempirical: depth | top1 (mass) | min-value node at depth | spine==min?")
spine = []
for d in range(15):
    L = layers[d]; tot = sum(L.values())
    (v1, c1), = L.most_common(1)
    spine.append(v1)
    mn = min_per_depth[d-2] if 2 <= d+1 <= 15 and d >= 1 else None
    mark = "==" if mn == v1 else ("!=" if mn else "")
    print(f"  {d+1:>2} | {v1:>6} ({c1/tot:.3f}) | {mn} {mark}")

print("\nspine value growth ratios:", [f"{b/a:.2f}" for a, b in zip(spine[1:], spine[2:])])

# --- corridor entropy
print("\nlayer entropy H_d (bits) and increments:")
prev = None
for d in range(15):
    L = layers[d]; tot = sum(L.values())
    H = -sum((c/tot)*math.log2(c/tot) for c in L.values())
    inc = f" (+{H-prev:.3f})" if prev is not None else ""
    prev = H
    if d % 2 == 0 or d > 10: print(f"  depth {d+1:>2}: H = {H:.3f}{inc}")
