"""R290-295: census on the bottom-up graph (79,545 nodes):
(1) TOWER LAW: in-degree of node v should equal sum over slide positions
    y=v*2^s of tern_ones(y) (stairs of every height 1..j land) -- verify the
    in-degree formula against v3-structure.
(2) TERNARY ALTERNATOR HUBS: are (9^j-1)/8 over-represented among top in-degree?
(3) CASTE LAYERING: in-degree by node mod 9 -- does the backward mint's
    caste system match the forward stationary law's mirror?"""
import sys
import networkx as nx
from collections import Counter
sys.stdout.reconfigure(encoding="utf-8")

G = nx.read_graphml("graph/collatz_backward.graphml")
# node ids are strings; value attribute
indeg = Counter()
val = {}
for n, d in G.nodes(data=True):
    val[n] = int(d["value"])
for n in G.nodes():
    indeg[val[n]] = G.in_degree(n)

alts = set((9**j - 1)//8 for j in range(1, 12))
print("(2) ternary alternators (9^j-1)/8 in the graph:")
ranked = sorted(indeg.items(), key=lambda t: -t[1])
rank_of = {v: i+1 for i, (v, _) in enumerate(ranked)}
for j in range(1, 8):
    a = (9**j - 1)//8
    if a in indeg:
        print(f"  j={j}: {a:>9,}  in-degree {indeg[a]:>3}  rank {rank_of[a]:>6,} / {len(ranked):,}")

import statistics
alt_deg = [indeg[a] for a in alts if a in indeg]
same_scale = []
for a in alts:
    if a in indeg:
        others = [indeg[v] for v in indeg if v != a and 0.5*a <= v <= 2*a]
        if others: same_scale.append((indeg[a], statistics.mean(others)))
print("  alternator in-degree vs same-scale mean:",
      [(d, f"{m:.2f}") for d, m in same_scale])

print("\n(3) mean in-degree by caste mod 9:")
by9 = Counter(); n9 = Counter()
for v, d in indeg.items():
    by9[v % 9] += d; n9[v % 9] += 1
for s in (1, 2, 4, 5, 7, 8, 0, 3, 6):
    if n9[s]: print(f"  mod9={s}: mean in-deg {by9[s]/n9[s]:.3f}  (n={n9[s]:,})")

print("\n(1) tower-law spot check: in-degree vs predicted stairs, 12 random nodes:")
import random
random.seed(290)
def tern_ones(n):
    j = 0
    while n % 3 == 1 and n > 0: j += 1; n //= 3
    return j
SLIDE_CAP = 3_000_000; VAL_CAP = 200_000
ok = bad = 0
for v in random.sample([x for x in indeg if x < 50_000], 200):
    pred = 0
    y, s = v, 0
    while y <= SLIDE_CAP:
        j = tern_ones(y)
        for k in range(1, j+1):
            Q = (y - (3**k - 1)//2) // 3**k
            n_ = (2*Q + 1) * 2**k - 1
            if n_ <= VAL_CAP and n_ != v: pred += 1
        y <<= 1; s += 1
    if pred == indeg[v]: ok += 1
    else: bad += 1
print(f"  exact matches {ok}/200, mismatches {bad} (mismatch = cap-truncation or dup)")
