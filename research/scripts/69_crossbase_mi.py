"""R661-670: CROSS-BASE MUTUAL INFORMATION along orbits. At T=0 CRT makes
n mod 3^b and w_1 = v2(3n+1) exactly independent. Along the orbit the map
mixes bases: measure I(n_T mod 9 ; w_{T+1}) exhaustively over n mod 2^12*3^5
(uniform start) for T = 0..8. Growth of MI = the base-coupling the Conversion
Thesis says carries the content."""
import sys, math
import numpy as np
from collections import Counter
sys.stdout.reconfigure(encoding="utf-8")

MOD = (1 << 12) * 3 ** 5
def v2(x):
    j = 0
    while x % 2 == 0: x //= 2; j += 1
    return j

for T in range(0, 9):
    joint = Counter(); m1 = Counter(); m2 = Counter(); tot = 0
    for n0 in range(1, MOD, 2):
        n = n0
        ok = True
        for _ in range(T):
            n = (3 * n + 1)
            n //= (1 << v2(n))
        if n % 3 == 0:  # springs: w defined anyway
            pass
        w = min(v2(3 * n + 1), 6)
        c = n % 9
        joint[(c, w)] += 1; m1[c] += 1; m2[w] += 1; tot += 1
    I = 0.0
    for (c, w), cnt in joint.items():
        p = cnt / tot
        I += p * math.log2(p * tot * tot / (m1[c] * m2[w]))
    print(f"T={T}: I(n_T mod 9 ; w_next) = {I:.6f} bits")
