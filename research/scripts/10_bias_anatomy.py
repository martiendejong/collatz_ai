"""E17: anatomy of the 2% borrow-chain bias — size dependence and digit-position localization.
Hypothesis: bias is a boundary effect of constrained low trits (3n+1 ends in trit 1, provably),
so it should localize at LSB positions and the whole-number average should decay ~ 1/log n."""
import sys, os, random, json
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "lib"))
sys.stdout.reconfigure(encoding="utf-8")
from collatz_lib import *
from collections import Counter

random.seed(31)
out = {}

def borrow_stats(x, pos_counter=None, tot_counter=None):
    """halve even x in base 3; count (r,r') transitions, optionally by digit distance from LSB."""
    ds = [int(c) for c in ternary(x)]
    L = len(ds)
    r = 0
    loc = Counter()
    for i, d in enumerate(ds):
        v = 3 * r + d
        r2 = v & 1
        loc[(r, r2)] += 1
        if pos_counter is not None:
            pos = L - 1 - i  # distance from least significant trit
            pos_counter[(min(pos, 30), r, r2)] += 1
        r = r2
    if tot_counter is not None:
        tot_counter.update(loc)
    return loc

# --- last-trit theorem check: x = 3n+1 always ends in trit 1 ---
bad = 0
for _ in range(20000):
    n = random.randrange(3, 1 << 40) | 1
    if (3 * n + 1) % 3 != 1: bad += 1
out["last_trit_of_3n+1_is_1"] = (bad == 0)

# --- bias vs size ---
size_bias = {}
for bits in (24, 32, 48, 64, 96, 128):
    c = Counter()
    for _ in range(1200):
        n = random.randrange(1 << (bits - 1), 1 << bits) | 1
        for _ in range(12):
            x = 3 * n + 1
            w = v2(x); y = x
            for _ in range(w):
                borrow_stats(y, tot_counter=c)
                y >>= 1
            n = y
            if n == 1: break
    p01 = c[(0, 1)] / (c[(0, 0)] + c[(0, 1)])
    p11 = c[(1, 1)] / (c[(1, 0)] + c[(1, 1)])
    size_bias[bits] = dict(P01=round(p01, 4), P11=round(p11, 4),
                           dev01=round(p01 - 1 / 3, 4), dev11=round(p11 - 2 / 3, 4))
out["bias_vs_size"] = size_bias

# --- bias by digit position (distance from LSB), 64-bit inputs ---
pos = Counter()
for _ in range(3000):
    n = random.randrange(1 << 63, 1 << 64) | 1
    for _ in range(12):
        x = 3 * n + 1
        w = v2(x); y = x
        for _ in range(w):
            borrow_stats(y, pos_counter=pos)
            y >>= 1
        n = y
        if n == 1: break
by_pos = {}
for p in list(range(0, 12)) + [20, 30]:
    d0 = pos[(p, 0, 0)] + pos[(p, 0, 1)]
    d1 = pos[(p, 1, 0)] + pos[(p, 1, 1)]
    if d0 > 500 and d1 > 500:
        by_pos[p] = dict(P01=round(pos[(p, 0, 1)] / d0, 4), P11=round(pos[(p, 1, 1)] / d1, 4))
out["bias_by_digit_position"] = by_pos

print(json.dumps(out, indent=1))
