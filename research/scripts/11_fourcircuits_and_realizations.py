"""E18: 4-circuit cycle exclusion (k_i <= 10, unbounded elements).
E19: integer-realization census of periodic index streams (2-adic positivity program)."""
import sys, os, json, math
sys.stdout.reconfigure(encoding="utf-8")
out = {}
LOG23 = math.log2(3)

def compose(shape):
    """closure constant for the cyclic chain; a1*(2^S - 3^K) = B."""
    A, B, S = 1, 0, 0
    m = len(shape)
    for i in range(m):
        k, w = shape[i]
        kn = shape[(i + 1) % m][0]
        A *= 3 ** k
        B = B * 3 ** k + ((1 << w) - 1) * (1 << S)
        S += w + kn
    return A, B, S

# ============ E18: 4-circuits, k_i <= 10, band 6 ============
KMAX, BAND = 10, 6
found4 = []
checked = 0
for k1 in range(1, KMAX + 1):
    for k2 in range(1, KMAX + 1):
        for k3 in range(1, KMAX + 1):
            for k4 in range(1, KMAX + 1):
                K = k1 + k2 + k3 + k4
                lo = max(math.ceil(K * LOG23), K + 4)
                for T in range(lo, lo + BAND):
                    W = T - K
                    if W < 4: continue
                    for w1 in range(1, W - 2):
                        for w2 in range(1, W - w1 - 1):
                            for w3 in range(1, W - w1 - w2):
                                w4 = W - w1 - w2 - w3
                                checked += 1
                                A, B, S = compose([(k1, w1), (k2, w2), (k3, w3), (k4, w4)])
                                D = (1 << S) - A
                                if D <= 0 or B % D: continue
                                a = B // D
                                if a > 0 and a % 2 == 1:
                                    found4.append((a, (k1, w1), (k2, w2), (k3, w3), (k4, w4)))
out["four_circuit"] = dict(checked=checked, kmax=KMAX,
    nontrivial=[f for f in found4 if f[0] != 1 or any(x != (1, 1) for x in f[1:])][:5],
    all_found=len(found4))

# ============ E19: integer realizations of periodic streams (pos AND neg) ============
# period-1: a = (2^w - 1)/(2^(k+w) - 3^k) -- any odd integer a (sign free) is a realized cycle.
real1 = []
for k in range(1, 41):
    for w in range(1, 41):
        D = (1 << (k + w)) - 3 ** k
        if D == 0: continue
        num = (1 << w) - 1
        if num % abs(D) == 0:
            a = num // D
            if a % 2 != 0 and a != 0:
                n = a * (1 << k) - 1
                real1.append((k, w, a, n))
out["period1_integer_realizations"] = real1

# period-2 realizations, k,w <= 18
real2 = []
for k1 in range(1, 19):
    for w1 in range(1, 19):
        for k2 in range(1, 19):
            for w2 in range(1, 19):
                A, B, S = compose([(k1, w1), (k2, w2)])
                D = (1 << S) - A
                if D == 0 or B == 0: continue
                if B % abs(D) == 0:
                    a = B // D
                    if a % 2 != 0 and a != 0:
                        n = a * (1 << k1) - 1
                        # skip pure repeats of period-1 patterns
                        if (k1, w1) == (k2, w2): continue
                        real2.append(((k1, w1), (k2, w2), a, n))
out["period2_integer_realizations_nonrepeat"] = real2[:20]
out["period2_count"] = len(real2)

print(json.dumps(out, indent=1))
