"""E22: period-4 realization census (k,w <= 7).
E23: climb-duration records — stopping-time champions and the E* tail exponent."""
import sys, os, json, math
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "lib"))
sys.stdout.reconfigure(encoding="utf-8")
from collatz_lib import v2
out = {}

def compose(shape):
    A, B, S = 1, 0, 0
    m = len(shape)
    for i in range(m):
        k, w = shape[i]
        kn = shape[(i + 1) % m][0]
        A *= 3 ** k
        B = B * 3 ** k + ((1 << w) - 1) * (1 << S)
        S += w + kn
    return A, B, S

# ============ E22: period-4 census ============
found = 0; nontriv = []
R = range(1, 8)
for k1 in R:
 for w1 in R:
  for k2 in R:
   for w2 in R:
    for k3 in R:
     for w3 in R:
      for k4 in R:
       for w4 in R:
        shape = ((k1, w1), (k2, w2), (k3, w3), (k4, w4))
        if shape[0] == shape[1] == shape[2] == shape[3]: continue
        A, B, S = compose(shape)
        D = (1 << S) - A
        if D == 0: continue
        if B % abs(D) == 0:
            a = B // D
            if a % 2 != 0 and a != 0:
                found += 1
                # exclude period-2 repeats ((x,y,x,y))
                if not (shape[0] == shape[2] and shape[1] == shape[3]):
                    nontriv.append((shape, a))
out["period4"] = dict(checked=7 ** 8, integer_realizations=found,
                      new_beyond_period2_repeats=nontriv[:8], new_count=len(nontriv))

# ============ E23: stopping-time records (T-map steps until value < start) ============
def stop_time(n, cap=1200):
    m = n; t = 0
    while t < cap:
        m = (3 * m + 1) // 2 if m % 2 else m // 2
        t += 1
        if m < n: return t
    return None  # censored

records = []
best = 0
NMAX = 5_000_000
for n in range(3, NMAX, 4):  # only n = 3 mod 4 can be records (k=1 drops fast)
    t = stop_time(n)
    if t is not None and t > best:
        best = t
        records.append((n, t))
out["stopping_time_records_n_t"] = records
out["known_literature_holders"] = [27, 703, 10087, 35655, 270271, 362343, 381727, 626331, 1027431, 4126387]

# tail exponent of stopping time distribution (E* observable)
from collections import Counter
import random
random.seed(51)
c = Counter()
M = 400000
for _ in range(M):
    n = random.randrange(1 << 30, 1 << 34) | 1
    t = stop_time(n, 500)
    c[t if t is not None else -1] += 1
# P(T >= t) tail slope in log2 between t=40 and t=120
def ge(t): return sum(v for key, v in c.items() if key == -1 or (key is not None and key >= t))
pts = [(t, math.log2(ge(t) / M)) for t in range(20, 141, 20) if ge(t) > 40]
slope = (pts[-1][1] - pts[0][1]) / (pts[-1][0] - pts[0][0])
out["tail_P(T>=t)_slope_bits_per_step"] = round(slope, 4)
out["needle_rate_F7"] = -0.05
print(json.dumps(out, indent=1))
