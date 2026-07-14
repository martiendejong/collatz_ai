"""R255: THE FUNNEL LAW. How does corridor concentration decay with height?
At depth d above 1 (in odd steps), measure the mass of the top-2 / top-2^?
most-visited odd values. The backward tree doubles per layer, springs kill 1/3;
question: does top-2 mass decay geometrically, and what is the corridor's
half-life?"""
import sys, random
from collections import Counter
sys.stdout.reconfigure(encoding="utf-8")
random.seed(255)

N = 150_000
DEPTH = 14
layers = [Counter() for _ in range(DEPTH)]
for _ in range(N):
    n = random.randrange(3, 1 << 44) | 1
    odds = [n]
    m = n
    while m != 1:
        m = 3*m+1 if m % 2 else m//2
        if m % 2 and m != 1: odds.append(m)
    for d in range(1, DEPTH+1):
        if len(odds) >= d:
            layers[d-1][odds[-d]] += 1

print(f"{'depth':>5} {'top1':>7} {'top2':>7} {'top4':>7} {'top8':>7} {'#distinct':>9}  top-2 values")
prev2 = None
for d in range(DEPTH):
    L = layers[d]; tot = sum(L.values())
    mc = L.most_common(8)
    t1 = mc[0][1]/tot
    t2 = sum(c for _, c in mc[:2])/tot
    t4 = sum(c for _, c in mc[:4])/tot
    t8 = sum(c for _, c in mc[:8])/tot
    ratio = f" (x{t2/prev2:.3f})" if prev2 else ""
    prev2 = t2
    print(f"{d+1:>5} {t1:>7.4f} {t2:>7.4f} {t4:>7.4f} {t8:>7.4f} {len(L):>9,}  {[v for v,_ in mc[:2]]}{ratio}")
