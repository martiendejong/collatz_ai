"""E12: are successive reloads independent along real orbits? (stochastic model test)
E13: 2-circuit cycle search (two-segment cycles)."""
import sys, os, random, json, math
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "lib"))
sys.stdout.reconfigure(encoding="utf-8")
from collatz_lib import *

random.seed(11)
out = {}

# ============ E12a: autocorrelation of successive (k_i) and (w_i) along orbits ============
ks, ws = [], []
for _ in range(4000):
    n = random.randrange(1 << 50, 1 << 60) | 1
    a, k = coords(n)
    for _ in range(40):
        a2, k2, w = macro_step(a, k)
        if a2 == 1: break
        ks.append((k, k2));
        a, k = a2, k2
def corr(pairs):
    xs = [p[0] for p in pairs]; ys = [p[1] for p in pairs]
    mx = sum(xs) / len(xs); my = sum(ys) / len(ys)
    cov = sum((x - mx) * (y - my) for x, y in pairs) / len(pairs)
    vx = sum((x - mx) ** 2 for x in xs) / len(xs)
    vy = sum((y - my) ** 2 for y in ys) / len(ys)
    return cov / math.sqrt(vx * vy)
out["corr_k_next_k"] = round(corr(ks), 5)
out["n_pairs"] = len(ks)

# conditional distribution: P(k'=1 | k) for k = 1..6 — should be 0.5 if memoryless
from collections import Counter
cond = {}
for k0 in range(1, 7):
    sub = [k2 for k, k2 in ks if k == k0]
    if len(sub) > 300:
        cond[k0] = round(sum(1 for x in sub if x == 1) / len(sub), 4)
out["P(k'=1 | k)"] = cond

# chi-square uniformity of a*3^k-1 mod 8 over orbit steps (drives reload law)
c = Counter()
for k, k2 in ks[:50000]:
    c[k2 if k2 < 5 else 5] += 1
tot = sum(c.values())
exp = {1: .5, 2: .25, 3: .125, 4: .0625, 5: .0625}
chi = sum((c[j] - exp[j] * tot) ** 2 / (exp[j] * tot) for j in exp)
out["chi2_reload_vs_geometric_df4"] = round(chi, 2)

# ============ E12b: pair (k_i, k_{i+2}) — longer-range dependence? ============
ks2 = []
for _ in range(2000):
    n = random.randrange(1 << 50, 1 << 60) | 1
    a, k = coords(n)
    hist = [k]
    for _ in range(40):
        a, k, w = macro_step(a, k)
        if a == 1: break
        hist.append(k)
    for i in range(len(hist) - 2):
        ks2.append((hist[i], hist[i + 2]))
out["corr_k_skip2"] = round(corr(ks2), 5)

# ============ E13: 2-circuit cycles ============
# chain: (a1,k1) -> (a2,k2) -> (a1,k1). Using exact transition:
#   m1 = (a1*3^k1 - 1)/2^w1 ; a2*2^k2 = m1+1
#   m2 = (a2*3^k2 - 1)/2^w2 ; a1*2^k1 = m2+1
# Search k1,k2,w1,w2 <= B, solve the linear system for a1 (rational), check integrality.
# Composing: a1*2^k1 = ((a2*3^k2 -1)/2^w2) + 1 and a2*2^k2 = ((a1*3^k1 -1)/2^w1) + 1
# => substitute: a1*(2^(k1+w1+w2) * 2^k2 ... do it with fractions.
from fractions import Fraction
B = 60
found = []
for k1 in range(1, B):
    for w1 in range(1, B):
        for k2 in range(1, B):
            for w2 in range(1, B):
                # a2 = (a1*3^k1 - 1 + 2^w1) / 2^(w1+k2)
                # a1 = (a2*3^k2 - 1 + 2^w2) / 2^(w2+k1)
                # linear: a1 = ((a1*3^k1 -1 + 2^w1)/2^(w1+k2) * 3^k2 - 1 + 2^w2) / 2^(w2+k1)
                # a1 * 2^(w2+k1) * 2^(w1+k2) = (a1*3^k1 -1 + 2^w1)*3^k2 + (2^w2 - 1)*2^(w1+k2)
                P = 1 << (w1 + w2 + k1 + k2)
                Q = 3 ** (k1 + k2)
                if P <= Q: continue
                num = (3 ** k2) * ((1 << w1) - 1) + ((1 << w2) - 1) * (1 << (w1 + k2))
                if num % (P - Q) == 0:
                    a1 = num // (P - Q)
                    if a1 > 0 and a1 % 2 == 1:
                        # verify it's a genuine 2-circuit (not the 1-circuit repeated)
                        a, k = coords(starter(a1, k1))
                        found.append((a1, k1, w1, k2, w2, starter(a1, k1)))
out["two_circuit_solutions_B60"] = found[:10]
out["two_circuit_count"] = len(found)

print(json.dumps(out, indent=1))
