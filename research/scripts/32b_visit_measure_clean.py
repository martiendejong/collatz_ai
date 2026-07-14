"""R258b: visit-measure law WITHOUT truncation bias: complete dyadic windows,
mean visit probability over ALL odd non-mult-3 values in window (zeros included),
plus the within-window distribution shape (spine vs bulk)."""
import sys, random, math
from collections import Counter
sys.stdout.reconfigure(encoding="utf-8")
random.seed(2582)

N = 120_000
CAP = 1 << 17
visits = Counter()
for _ in range(N):
    n = random.randrange(3, 1 << 40) | 1
    m = n
    while m != 1:
        m = 3*m+1 if m % 2 else m//2
        if m % 2 and m < CAP: visits[m] += 1

print("complete dyadic windows (all odd v, not mult 3, zeros included):")
print(f"{'window':>18} {'mean P':>10} {'exponent step':>13} {'max/mean':>9} {'P(0 visits)':>11}")
prev = None
for e in range(7, 17):
    lo, hi = 1 << e, 1 << (e+1)
    vals = [v for v in range(lo | 1, hi, 2) if v % 3 != 0]
    ps = [visits.get(v, 0)/N for v in vals]
    mp = sum(ps)/len(ps)
    step = f"{math.log(mp/prev)/math.log(2):+.3f}" if prev else "     -"
    prev = mp
    zfrac = sum(1 for p in ps if p == 0)/len(ps)
    print(f"[2^{e:>2}, 2^{e+1:>2}) {mp:>10.2e} {step:>13} {max(ps)/mp:>9.1f} {zfrac:>11.3f}")
print("\nexponent step = log2(mean P ratio) per doubling; -1 predicted by uniform-window law")
