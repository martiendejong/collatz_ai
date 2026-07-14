"""R421-432: CYCLE CENSUS through period 12 (odd steps k <= 12), positive AND
negative integers. Cycle equation: n0*(2^s - 3^k) = sum_i 3^{k-1-i} * 2^{S_i},
S_i = w_1+...+w_i (S_0=0), over all compositions (w_1..w_k), w_i>=1.
Positive cycles need 2^s > 3^k; negative need 2^s < 3^k. Known: {1} positive;
{-1,-5,-17} negative. Verify completeness through k=12."""
import sys, math
from itertools import combinations
sys.stdout.reconfigure(encoding="utf-8")

found = set()
for k in range(1, 13):
    smin, smax = k, int(k * math.log2(3)) + 2
    for s in range(smin, smax + 1):
        D = (1 << s) - 3 ** k
        if D == 0: continue
        # compositions of s into k parts >= 1 <-> subsets: S_1<...<S_{k-1} strictly in [1,s-1]
        for cuts in combinations(range(1, s), k - 1):
            S = (0,) + cuts
            num = sum(3 ** (k - 1 - i) * (1 << S[i]) for i in range(k))
            if num % D == 0:
                n0 = num // D
                if n0 % 2 == 1 and n0 != 0:
                    # canonical: min |odd| element of the cycle
                    orb = [n0]; m = n0
                    for w in [S[i+1]-S[i] for i in range(k-1)] + [s - S[k-1]]:
                        m = (3 * m + 1) >> w
                        orb.append(m)
                    found.add(min(orb[:-1], key=lambda x: (abs(x), x)))
print("all integer Syracuse cycles with <=12 odd steps (canonical min element):")
print(sorted(found, key=lambda x: (abs(x), x)))
print("\nknown: 1 (positive), -1, -5, -17 (negative). New ones would appear above.")
