"""E1: closed-form (F,P) transition, the exact rewriting theorem, reload structure."""
import sys, os, random, json
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "lib"))
sys.stdout.reconfigure(encoding="utf-8")
from collatz_lib import *
from collections import Counter

random.seed(1)
out = {}

# --- 1. closed-form macro transition verified against raw iteration ---
bad = 0
for _ in range(50000):
    n = random.randrange(3, 1 << 48) | 1
    a, k = coords(n)
    a2, k2, w = macro_step(a, k)
    # raw: run k steps of (3n+1)/2 then strip 2s
    m = n
    for _ in range(k):
        m = (3 * m + 1) // 2
    assert m == a * 3 ** k - 1
    m >>= v2(m)
    if coords(m) != (a2, k2):
        bad += 1
out["transition_verified"] = bad == 0

# --- 2. THE REWRITING THEOREM: binary word a|1^k  ->  ternary word (a-1)|2^k ---
# starter = a*2^k - 1 : binary = bits of a (minus its final 1 merging into the run? no: exact claim below)
# endpoint = a*3^k - 1 = (a-1)*3^k + (3^k - 1): ternary = digits of (a-1) followed by k twos.
ok = True
examples = []
for a in range(1, 400, 2):
    for k in range(1, 12):
        endpoint = a * 3 ** k - 1
        want = (ternary(a - 1) if a > 1 else "") + "2" * k
        got = ternary(endpoint)
        if got != want:
            ok = False
        if len(examples) < 6 and a in (1, 5, 7) and k in (1, 3):
            examples.append((starter(a, k), bin(starter(a, k))[2:], endpoint, got))
out["rewriting_theorem"] = ok
out["rewriting_examples"] = examples

# --- 3. reload (w and next k) distributions: geometric? ---
wc, kc, tot = Counter(), Counter(), 0
for _ in range(300000):
    n = random.randrange(1 << 30, 1 << 44) | 1
    a, k = coords(n)
    a2, k2, w = macro_step(a, k)
    wc[min(w, 12)] += 1
    kc[min(k2, 12)] += 1
    tot += 1
out["w_dist"] = {j: round(wc[j] / tot, 5) for j in range(1, 9)}
out["k2_dist"] = {j: round(kc[j] / tot, 5) for j in range(1, 9)}

# --- 4. reload predictability depth: v2(a*3^k-1) determined by (a mod 2^j, k mod 2^(j-2)) ---
depth_ok = {}
for j in (4, 6, 8, 10):
    good = True
    for _ in range(4000):
        a1 = random.randrange(1, 1 << 22) * 2 + 1
        k1 = random.randrange(1, 80)
        a2_ = a1 + (random.randrange(1, 500) << j)
        k2_ = k1 + (1 << max(j - 2, 1)) * random.randrange(1, 6)
        if min(v2(a1 * 3 ** k1 - 1), j) != min(v2(a2_ * 3 ** k2_ - 1), j):
            good = False
            break
    depth_ok[j] = good
out["reload_predictability"] = depth_ok

print(json.dumps(out, indent=1))
