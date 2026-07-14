"""R272-277: MULTIFRACTAL SPECTRUM of the visit measure.
Within each dyadic window [2^e, 2^{e+1}), normalize the visit distribution and
compute generalized participation exponents:
  D_q = (1/(1-q)) * log(sum p_i^q) / log(N_window)
D_1 (information dim, via entropy), D_2, D_4, D_inf (max). If the measure were
uniform in-window: all D_q = 1. Multifractality = spread of D_q, and its trend
in e tells whether concentration sharpens with scale."""
import sys, random, math
from collections import Counter
sys.stdout.reconfigure(encoding="utf-8")
random.seed(272)

N = 250_000
CAP = 1 << 18
visits = Counter()
for _ in range(N):
    n = random.randrange(3, 1 << 44) | 1
    m = n
    while m != 1:
        m = 3*m+1 if m % 2 else m//2
        if m % 2 and m < CAP: visits[m] += 1

print(f"{'window':>14} {'Nw':>6} {'D1':>6} {'D2':>6} {'D4':>6} {'Dinf':>6} {'top-1 share':>11}")
for e in range(8, 18):
    lo, hi = 1 << e, 1 << (e+1)
    vals = [v for v in range(lo | 1, hi, 2) if v % 3 != 0]
    cs = [visits.get(v, 0) for v in vals]
    tot = sum(cs)
    if tot < 5000: continue
    ps = [c/tot for c in cs if c > 0]
    Nw = len(vals); lnN = math.log(Nw)
    H = -sum(p*math.log(p) for p in ps)
    D1 = H/lnN
    D2 = -math.log(sum(p*p for p in ps))/lnN
    D4 = -math.log(sum(p**4 for p in ps))/(3*lnN)
    Di = -math.log(max(ps))/lnN
    print(f"[2^{e:>2},2^{e+1:>2}) {Nw:>6} {D1:>6.3f} {D2:>6.3f} {D4:>6.3f} {Di:>6.3f} {max(ps):>11.4f}")
print("\nuniform measure would give D_q = 1 for all q; spread = multifractality")
