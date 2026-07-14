"""E20: the Sign Theorem + Catalan structure of the known cycles.
E21: period-3 realization census (should find nothing new -- there are only 3 cycles on Z)."""
import sys, os, json, math
sys.stdout.reconfigure(encoding="utf-8")
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

# ============ E20a: verify B > 0 always and sign(a1) = sign(2^S - 3^K) ============
import random
random.seed(41)
ok = True
for _ in range(20000):
    m = random.randint(1, 5)
    shape = [(random.randint(1, 12), random.randint(1, 12)) for _ in range(m)]
    A, B, S = compose(shape)
    if B <= 0: ok = False; break
out["B_always_positive"] = ok
out["sign_theorem"] = "a1 = B/(2^S - 3^K), B>0  =>  sign(a1) = sign(2^(K+W) - 3^K): positive cycles need net-falling shapes (W > 0.585 K)"

# ============ E20b: the unit-gap (Catalan) table: min |2^S - 3^K| for K <= 40 ============
gaps = []
for K in range(1, 26):
    best = None
    for S in (math.floor(K * math.log2(3)), math.ceil(K * math.log2(3))):
        g = (1 << S) - 3 ** K
        if best is None or abs(g) < abs(best[1]):
            best = (S, g)
    gaps.append((K, best[0], best[1]))
out["min_gap_table_K_S_gap"] = gaps[:14]
unit = [(K, S, g) for K, S, g in gaps if abs(g) == 1]
out["unit_gaps"] = unit
out["catalan_note"] = "unit gaps only at (K,S)=(1,2): 4-3=+1 -> trivial cycle; (2,3): 8-9=-1 -> the -5 cycle. Mihailescu: no others ever."

# known cycles vs their gaps:
out["known_cycles_gap"] = {
    "n=1  shape (1,1)":  "2^2 - 3^1 = +1  (a = 1/1 = 1)",
    "n=-5 shape (2,1)":  "2^3 - 3^2 = -1  (a = 1/-1 = -1)",
    "n=-17 shape (4,1)(3,3)": f"2^11 - 3^7 = {2**11 - 3**7}  (divisibility: B=|gap|*|a|)",
}

# ============ E21: period-3 census, k,w <= 10 ============
found = []
for k1 in range(1, 11):
 for w1 in range(1, 11):
  for k2 in range(1, 11):
   for w2 in range(1, 11):
    for k3 in range(1, 11):
     for w3 in range(1, 11):
        shape = [(k1, w1), (k2, w2), (k3, w3)]
        A, B, S = compose(shape)
        D = (1 << S) - A
        if D == 0: continue
        if B % abs(D) == 0:
            a = B // D
            if a % 2 != 0 and a != 0:
                # exclude repeats of period-1 (all three equal)
                if (k1, w1) == (k2, w2) == (k3, w3): continue
                found.append((shape, a, a * (1 << k1) - 1))
out["period3_nonrepeat_realizations"] = found[:12]
out["period3_count"] = len(found)
print(json.dumps(out, indent=1))
