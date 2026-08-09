# 352: the price of "eventually" in the backward tree.
# Martien: "we only still need to determine that every number eventually appears
# in the backward tree" — which is verbatim the conjecture. The backward tree
# from 1 reaches every verified number, but through DETOURS (27 needs 9232).
# Quantify: F(b) = P(peak(n)/n >= 2^b) over random odd n.
# Prediction from the growth law (Obs 534/538): the cheapest way to gain bits is
# a k-ones run (prob 2^-k, gain k*(log2(3)-1) bits), so the large-deviation rate
# is theta = 1/(log2(3)-1) = 1.7095:  F(b) ~ C * 2^(-theta*b).
import random
from math import log2

random.seed(352)
N_SAMPLES = 120000
LO, HI = 10**5, 10**7

ratios = []
for _ in range(N_SAMPLES):
    n = random.randrange(LO | 1, HI, 2)
    m, peak = n, n
    while m != 1:
        m = 3 * m + 1 if m % 2 else m // 2
        if m > peak:
            peak = m
    ratios.append(log2(peak) - log2(n))

ratios.sort()
print(f"{N_SAMPLES} random odd n in [{LO},{HI}]: excursion ratio log2(peak/n)")
print(f"  median {ratios[len(ratios)//2]:.2f} bits, max {ratios[-1]:.2f} bits")
print()
print("  b   P(ratio>=b)      log2 P    slope (local theta)")
prev = None
import bisect
for b in range(2, 22, 2):
    cnt = len(ratios) - bisect.bisect_left(ratios, b)
    if cnt < 10:
        break
    p = cnt / len(ratios)
    l2 = log2(p)
    slope = (prev - l2) / 2 if prev is not None else float('nan')
    print(f"  {b:>2}  {p:11.6f}  {l2:9.3f}   {slope:6.3f}")
    prev = l2
print()
print(f"prediction: theta = 1/(log2(3)-1) = {1/(log2(3)-1):.4f}")
print("(local slopes should approach this from below as b grows)")
