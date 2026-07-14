"""R258: THE VISIT MEASURE LAW. Why does the funnel spine plateau (R255)?
Hypothesis: P(random orbit visits odd v) ~ c/v (basin density scale-invariance).
Then the top node at each depth is simply the smallest node reachable at that
depth, and its mass doesn't dilute. Test: measure visit frequency of many odd
v across random orbits, regress log P against log v."""
import sys, random, math
from collections import Counter
sys.stdout.reconfigure(encoding="utf-8")
random.seed(258)

N = 120_000
visits = Counter()
for _ in range(N):
    n = random.randrange(3, 1 << 40) | 1
    m = n
    while m != 1:
        m = 3*m+1 if m % 2 else m//2
        if m % 2 and m < 100_000: visits[m] += 1

# regress log(visit P) on log v for well-sampled v
pts = [(math.log(v), math.log(c/N)) for v, c in visits.items() if c >= 50 and v % 3 != 0]
n_ = len(pts)
sx = sum(x for x,_ in pts); sy = sum(y for _,y in pts)
sxx = sum(x*x for x,_ in pts); sxy = sum(x*y for x,y in pts)
slope = (n_*sxy - sx*sy)/(n_*sxx - sx*sx)
inter = (sy - slope*sx)/n_
print(f"visit-measure regression over {n_} odd values (v<1e5, v not mult of 3):")
print(f"  P(visit v) ~ {math.exp(inter):.3f} * v^{slope:.4f}   (hypothesis: exponent -1)")
# residual structure by caste
for r in (1, 2):
    pr = [(math.log(v), math.log(c/N)) for v, c in visits.items() if c >= 50 and v % 3 == r]
    m_ = len(pr)
    sx2 = sum(x for x,_ in pr); sy2 = sum(y for _,y in pr)
    sxx2 = sum(x*x for x,_ in pr); sxy2 = sum(x*y for x,y in pr)
    sl = (m_*sxy2 - sx2*sy2)/(m_*sxx2 - sx2*sx2)
    print(f"  caste v=={r} mod 3: exponent {sl:.4f} ({m_} values)")
# sample devations: most over/under-visited relative to fit
def resid(v, c): return math.log(c/N) - (inter + slope*math.log(v))
top = sorted(((resid(v,c), v) for v,c in visits.items() if c >= 200), reverse=True)[:5]
bot = sorted(((resid(v,c), v) for v,c in visits.items() if c >= 200))[:5]
print("  most OVER-visited vs law:", [(v, f"{r:+.2f}") for r, v in top])
print("  most UNDER-visited vs law:", [(v, f"{r:+.2f}") for r, v in bot])
