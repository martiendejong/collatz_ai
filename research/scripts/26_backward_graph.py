"""R219-222: build the bottom-up Collatz graph via the Backward Stairway Rule.
From v=1: slide upward (doublings); at each slide position, read the TERNARY trailing ones j
of y = m*2^s; for each k<=j a staircase of height k lands here, from n = (2Q+1)*2^k - 1,
Q = (y - (3^k-1)/2)/3^k. Nodes = odd numbers; edge n->m = one macro-step (attrs k, w).
Outputs: GraphML + GEXF + stats + a rendered PNG."""
import sys, math
import networkx as nx
sys.stdout.reconfigure(encoding="utf-8")

VAL_CAP = 200_000        # only keep nodes with value <= cap
SLIDE_CAP = 3_000_000    # explore slide positions y = m*2^s up to this
OUT = "graph/collatz_backward"

def tern(n):
    s = ""
    while n: s = str(n % 3) + s; n //= 3
    return s or "0"
def tern_ones(n):
    j = 0
    while n % 3 == 1 and n > 0: j += 1; n //= 3
    return j
def trailing_ones_bin(n):
    k = 0
    while n & 1: n >>= 1; k += 1
    return k

G = nx.DiGraph()
def node_attrs(n, depth):
    k = trailing_ones_bin(n)
    a = (n + 1) >> k
    k0 = 1 if a % 4 == 1 else 2
    y = (a * 3 ** k0 - 1) // 2
    M = 3 * y + 1
    return dict(value=n, binary=bin(n)[2:], ternary=tern(n), family=int(a), seqindex=int(k),
                ftype=int(0 if M % 4 == 0 else 2), depth=int(depth))

from collections import deque
G.add_node(1, **node_attrs(1, 0))
q = deque([(1, 0)])
seen = {1}
edges = 0
while q:
    m, d = q.popleft()
    s = 0
    y = m
    while y <= SLIDE_CAP:
        j = tern_ones(y)
        for k in range(1, j + 1):
            Q = (y - (3 ** k - 1) // 2) // 3 ** k
            n = (2 * Q + 1) * 2 ** k - 1
            if n <= VAL_CAP and n != m and n not in seen:
                seen.add(n)
                G.add_node(n, **node_attrs(n, d + 1))
                G.add_edge(n, m, k=int(k), w=int(s + 1))
                edges += 1
                q.append((n, d + 1))
            elif n <= VAL_CAP and n in seen and not G.has_edge(n, m):
                pass  # each n has exactly one forward edge; first discovery is the true one
        y <<= 1; s += 1

print(f"graph built: {G.number_of_nodes():,} nodes, {G.number_of_edges():,} edges (a TREE rooted at 1)")

import os
os.makedirs("graph", exist_ok=True)
nx.write_graphml(G, OUT + ".graphml")
nx.write_gexf(G, OUT + ".gexf")
print(f"written: {OUT}.graphml and .gexf")

# ===== what do we SEE =====
indeg = dict(G.in_degree())
hubs = sorted(indeg.items(), key=lambda t: -t[1])[:10]
print("\ntop landing-hubs (in-degree = stairs landing on this node's slide within caps):")
for v, dg in hubs:
    print(f"  {v:>7}: in-degree {dg:>3}   ternary {tern(v)}")
from collections import Counter
dc = Counter(indeg.values())
print("\nin-degree distribution:", dict(sorted(dc.items())[:10]))
depths = [d for _, d in G.nodes(data="depth")]
print(f"depth (macro-steps to 1): max {max(depths)}, mean {sum(depths)/len(depths):.1f}")
kdist = Counter(k for _, _, k in G.edges(data="k"))
tot = sum(kdist.values())
print("stair-height distribution on edges:", {k: round(v/tot, 4) for k, v in sorted(kdist.items())[:6]})

# ===== render =====
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sub = G.subgraph([n for n, d in G.nodes(data="depth") if d <= 14])
pos = {}
import collections
bydepth = collections.defaultdict(list)
for n, d in sub.nodes(data="depth"): bydepth[d].append(n)
for d, ns in bydepth.items():
    ns.sort()
    for i, n in enumerate(ns):
        pos[n] = (i - len(ns)/2, -d)
plt.figure(figsize=(14, 9))
cols = ["#1a7a4a" if t == 0 else "#8b1e3f" for _, t in sub.nodes(data="ftype")]
nx.draw_networkx_edges(sub, pos, alpha=0.25, arrows=False, width=0.5)
nx.draw_networkx_nodes(sub, pos, node_size=8, node_color=cols)
plt.title(f"Collatz backward graph, built bottom-up by the Ternary-Tail Rule ({sub.number_of_nodes()} nodes, depth<=14)\n"
          "green = type-0 families (x9+1 chain), red = type-2 | y = macro-steps to 1")
plt.axis("off"); plt.tight_layout()
plt.savefig("figures/backward_graph.png", dpi=130)
print("figure: figures/backward_graph.png")
