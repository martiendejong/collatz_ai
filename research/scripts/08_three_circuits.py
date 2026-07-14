"""E14: 3-circuit cycle exclusion via pruned exact search.
Chain: (a1,k1)->(a2,k2)->(a3,k3)->(a1,k1). Composing the exact transition:
  a1*(2^(K+W) - 3^K) = B_total,  K=k1+k2+k3, W=w1+w2+w3
where B_total accumulates iteratively. Prune: a1>=1 odd requires 2^(K+W) barely above 3^K.
"""
import sys, os, json, math
sys.stdout.reconfigure(encoding="utf-8")
out = {}

LOG23 = math.log2(3)
KMAX = 24          # per-segment k_i <= KMAX
BAND = 8           # (K+W) - ceil(K log2 3) within band
found = []
checked = 0

pow3 = [3 ** i for i in range(3 * KMAX + 1)]

def compose(shape):
    """shape = [(k1,w1),(k2,w2),(k3,w3)] with k_{i+1} consumed after w_i.
    a_{i+1} = (a_i*3^{k_i} - 1 + 2^{w_i}) / 2^{w_i + k_{i+1}}
    Maintain a_{i} = (A*a1 + B)/2^S."""
    A, B, S = 1, 0, 0
    m = len(shape)
    for i in range(m):
        k, w = shape[i]
        kn = shape[(i + 1) % m][0]
        # a' = (a*3^k - 1 + 2^w)/2^(w+kn); with a = (A a1 + B)/2^S:
        # a' = (A*3^k a1 + B*3^k + (2^w - 1)*2^S) / 2^(S + w + kn)
        A = A * pow3[k]
        B = B * pow3[k] + ((1 << w) - 1) * (1 << S)
        S = S + w + kn
    return A, B, S   # closure: a1 = (A a1 + B)/2^S -> a1*(2^S - A) = B

for k1 in range(1, KMAX + 1):
    for k2 in range(1, KMAX + 1):
        for k3 in range(1, KMAX + 1):
            K = k1 + k2 + k3
            lo = math.ceil(K * LOG23)
            for T in range(max(lo, K + 3), lo + BAND):  # T = K + W, W = T-K >= 3
                W = T - K
                if W < 3: continue
                # compositions W = w1+w2+w3, wi>=1
                for w1 in range(1, W - 1):
                    for w2 in range(1, W - w1):
                        w3 = W - w1 - w2
                        checked += 1
                        A, B, S = compose([(k1, w1), (k2, w2), (k3, w3)])
                        D = (1 << S) - A
                        if D <= 0 or B % D: continue
                        a1 = B // D
                        if a1 > 0 and a1 % 2 == 1:
                            found.append((a1, (k1, w1), (k2, w2), (k3, w3)))

out["shapes_checked"] = checked
out["k_max"] = KMAX
out["band"] = BAND
# filter: trivial = all segments (1,1) with a=1
nontrivial = [f for f in found if f[0] != 1 or any(kw != (1, 1) for kw in f[1:])]
out["solutions"] = found[:10]
out["nontrivial_solutions"] = nontrivial[:10]
print(json.dumps(out, indent=1))
