# 401: the FLOW MAP of the funnel — fraction of random orbits passing each
# small odd station; the graph structure is w -> syr(w) with flow weights.
# Also test the harmonic hypothesis P(pass w) ~ C/w, and export JSON for the
# course visualization.
import random, json
from collections import Counter
random.seed(401)

def syr(n):
    m = 3*n + 1
    return m >> ((m & -m).bit_length() - 1)

N = 6000
LIM = 1 << 22
cnt = Counter()
for _ in range(N):
    n = random.getrandbits(64) | 1
    seen = set()
    while n != 1:
        if n < LIM and n not in seen:
            seen.add(n)
            cnt[n] += 1
        n = syr(n)
    cnt[1] += 1

flow = {w: c/N for w, c in cnt.items()}
top = sorted(flow.items(), key=lambda kv: -kv[1])[:36]
print("top-stations (P(baan passeert w)):")
for w, f in top[:20]:
    print(f"  w={w:>6}: {f:.4f}")

# harmonic test: P ~ C/w over deciles of w
import math
ws = sorted(w for w in flow if flow[w] > 0.01 and w > 4)
print("\nharmonische toets (P*w moet ~constant zijn op de hoofdtak):")
for w in [5, 11, 17, 53, 161, 485, 1457]:
    if w in flow:
        print(f"  w={w:>5}: P = {flow[w]:.4f}  P*w = {flow[w]*w:.2f}")

# export graph data: nodes with flow >= 0.02, edges w -> syr(w)
nodes = {w: f for w, f in flow.items() if f >= 0.02}
# distance-to-1 for layout
def dist1(w):
    d = 0
    while w != 1:
        w = syr(w); d += 1
    return d
data = []
for w, f in sorted(nodes.items()):
    data.append({"w": w, "f": round(f, 4), "to": syr(w) if w != 1 else 0, "d": dist1(w)})
out = json.dumps(data)
open(r"E:\projects\collatz\research\flow_map.json", "w").write(out)
print(f"\n{len(data)} knopen met flow >= 2% geexporteerd; maximale afstand tot 1: "
      f"{max(x['d'] for x in data)}")
