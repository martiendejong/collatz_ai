# 383: the pattern of the FRONT bits under pure repeated x3 (Martien's question).
# Claims to verify:
#  (H1) bit-length increments of 3^k form the STURMIAN word of slope
#       alpha-1 = 0.585: d_k = 2 iff {k*alpha} >= 2-alpha. Exact match?
#  (H2) leading-bit blocks follow the base-2 Benford law:
#       P(front "10") = log2(3/2) = 0.585, P("11") = 2-alpha = 0.415;
#       3-bit blocks 100/101/110/111: log2(5/4), log2(6/5), log2(7/6), log2(8/7).
#  (H3) the front is a ROTATION: leading j bits of 3^k depend only on
#       frac(k*alpha) to precision 2^-j (zero entropy, fully predictable).
from math import log2, floor
import numpy as np

ALPHA = log2(3.0)
N = 20000
x = 1
lens = []
lead4 = []
fracs = []
for k in range(1, N+1):
    x *= 3
    L = x.bit_length()
    if L < 4:
        continue
    lens.append(L)
    lead4.append(x >> (L-4))       # top-4-bit value in [8,16)
    fracs.append((k*ALPHA) % 1.0)
lens = np.array(lens); lead4 = np.array(lead4); fracs = np.array(fracs)

# H1: increments vs Sturmian rule
d = np.diff(lens)
rule = (fracs[:-1] >= 2 - ALPHA).astype(int) + 1
match = (d == rule).mean()
print(f"(H1) lengte-incrementen: {np.bincount(d)[1:]} (1'en/2'en); "
      f"Sturmian-regel d=2 <=> {{k*alpha}} >= {2-ALPHA:.4f}: match {match:.6f}")
print(f"     fractie 2'en = {(d == 2).mean():.5f}  (alpha-1 = {ALPHA-1:.5f})")

# H2: Benford-2 for 2- and 3-bit fronts
p10 = ((lead4 >= 8) & (lead4 < 12)).mean()
p11 = (lead4 >= 12).mean()
print(f"(H2) P(front 10) = {p10:.5f} (Benford: {log2(1.5):.5f})   "
      f"P(front 11) = {p11:.5f} (Benford: {2-ALPHA:.5f})")
for v in range(8, 16):
    p = (lead4 == v).mean()
    b = log2((v+1)/v)
    print(f"     front {v:04b}: {p:.5f}  (Benford log2({v+1}/{v}) = {b:.5f})")

# H3: rotation property — same frac(k*alpha) bucket => same leading bits
buckets = {}
coll = 0; tot = 0
for k in range(len(fracs)):
    b = int(fracs[k]*4096)
    if b in buckets:
        tot += 1
        if buckets[b] != lead4[k]:
            coll += 1
    else:
        buckets[b] = lead4[k]
print(f"(H3) rotatie-eigenschap: {tot} bucket-herhalingen (precisie 2^-12), "
      f"{coll} met andere top-4-bits ({100*coll/max(tot,1):.2f}% · verwacht ~0 op randen na)")
