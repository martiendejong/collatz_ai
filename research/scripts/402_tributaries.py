# 402: TRIBUTARIES — inflow decomposition of the funnel stations and the flow
# profile ALONG the artery upward (extend the 27-finale chain by max-flow
# predecessors; test the decay law). Edge flows counted from fresh sample.
import random
from collections import Counter
random.seed(402)

def syr(n):
    m = 3*n + 1
    return m >> ((m & -m).bit_length() - 1)

N = 6000
LIM = 1 << 24
edge = Counter()
node = Counter()
for _ in range(N):
    n = random.getrandbits(64) | 1
    seen = set()
    prev = None
    while n != 1:
        if n < LIM and n not in seen:
            seen.add(n); node[n] += 1
            if prev is not None:
                edge[(prev, n)] += 1
            prev = n
        else:
            prev = n if n < LIM else None
        n = syr(n)
    node[1] += 1
    if prev is not None:
        edge[(prev, 1)] += 1

def inflow(w, top=4):
    ins = [(m, c) for (m, t), c in edge.items() if t == w]
    ins.sort(key=lambda x: -x[1])
    tot = sum(c for _, c in ins)
    return tot, ins[:top]

print("instroom-decompositie van de trechterstations:")
for w in [1, 5, 13, 53, 11, 17]:
    tot, ins = inflow(w)
    parts = ", ".join(f"{m}:{c/N:.3f}" for m, c in ins)
    print(f"  w={w:>4} (totaal {tot/N:.3f}): {parts}")

# de ader omhoog: volg vanaf 5 telkens de voorganger met maximale flow
print("\nde ader omhoog (max-flow-voorganger vanaf 5):")
w = 5
chain = [(5, node[5]/N)]
for _ in range(22):
    tot, ins = inflow(w, top=1)
    if not ins or ins[0][1] < 30:
        break
    w = ins[0][0]
    chain.append((w, node[w]/N))
import math
print("  " + " <- ".join(f"{w}({f:.3f})" for w, f in chain))
# decay fit along artery
ws = [w for w, f in chain if w > 5]
fs = [f for w, f in chain if w > 5]
if len(ws) > 4:
    lx = [math.log(w) for w in ws]; ly = [math.log(f) for f in fs]
    n_ = len(lx)
    sx = sum(lx); sy = sum(ly); sxx = sum(x*x for x in lx); sxy = sum(x*y for x, y in zip(lx, ly))
    slope = (n_*sxy - sx*sy)/(n_*sxx - sx*sx)
    print(f"\n  machtwet langs de ader: flow ~ w^{slope:.3f}")
