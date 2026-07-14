"""E31: class-count K-L system (correct semantics).
pi_(b,a,j)(x) = count of P-elements <= x with parity b, = a mod 3^j.
Identities: (0,a,j): pi = sum_{b'} pi_(b',a/2,j)(x/2)         [exact]
            (1,a,j): pi = pi_(0,3a+1,j+1)(3x+1)               [exact, deeper modulus]
Unfold tree with optimal stopping (harvest when net scale 3^o/2^e < 1), leaves are class
variables at their own level. Extra operator term: u_coarse >= sum of refinements (consistency).
Certificate: u > 0 with u <= Op(u) => pi_c(x) >= eps*u_c*x^gamma for all classes, by induction.
"""
import sys, json, math
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
out = {}
L23 = math.log2(3)

def solve(K, T, piters=40):
    # pack states (b, a mod 3^j, j) for j=0..K
    offs = []; tot = 0
    for j in range(K + 1):
        offs.append(tot); tot += 2 * 3 ** j
    def idx(b, a, j): return offs[j] + b * 3 ** j + a
    # transitions
    evenA = np.full(tot, -1, dtype=np.int64)   # (0,a,j) -> (b', a*inv2, j) both parities
    oddA = np.full(tot, -1, dtype=np.int64)    # (1,a,j) -> (0, 3a+1, j+1)
    j_of = np.empty(tot, dtype=np.int64); b_of = np.empty(tot, dtype=np.int64)
    for j in range(K + 1):
        M = 3 ** j
        inv2 = pow(2, -1, M) if M > 1 else 0
        for b in (0, 1):
            for a in range(M):
                i = idx(b, a, j); j_of[i] = j; b_of[i] = b
                if b == 0:
                    evenA[i] = idx(0, (a * inv2) % M if M > 1 else 0, j)  # parity-0 child; parity-1 = +M
                else:
                    if j < K:
                        oddA[i] = idx(0, (3 * a + 1) % (3 * M), j + 1)
    def apply_op(gamma, u):
        w1 = 2.0 ** (-gamma); w2 = 3.0 ** gamma
        # V-DP over (state, e, o): value normalized by s^gamma (s = 3^o/2^e)
        # W(state,e,o) = max(harvest: u[state] if 3^o < 2^e, even: w1*(W(c0,e+1,o)+W(c1,e+1,o)),
        #                    odd: w2*W(child,e,o+1))
        W = np.zeros((tot, T + 2, K + 2))
        for o in range(K, -1, -1):
            for e in range(T, -1, -1):
                col = np.zeros(tot)
                harv_ok = (3.0 ** o) < (2.0 ** e)
                if harv_ok:
                    col = u.copy()
                if e < T:
                    ev = evenA >= 0
                    cont = np.zeros(tot)
                    cont[ev] = w1 * (W[evenA[ev], e + 1, o] + W[evenA[ev] + 3 ** j_of[ev] * 0 + (3 ** j_of[ev]), e + 1, o])
                    col = np.maximum(col, cont)
                if o < K:
                    od = oddA >= 0
                    cont2 = np.zeros(tot)
                    cont2[od] = w2 * W[oddA[od], e, o + 1]
                    col = np.maximum(col, cont2)
                W[:, e, o] = col
        V = W[:, 0, 0]
        # consistency term: u_coarse >= sum over refinements at j+1 (same parity)
        newu = V.copy()
        for j in range(K - 1, -1, -1):
            M = 3 ** j
            for b in (0, 1):
                fine = newu[offs[j + 1] + b * 3 * M: offs[j + 1] + b * 3 * M + 3 * M].reshape(3, M)
                seg = slice(offs[j] + b * M, offs[j] + (b + 1) * M)
                newu[seg] = np.maximum(newu[seg], fine.sum(axis=0))
        return newu
    def rho(gamma):
        u = np.ones(tot)
        r = 0.0
        for it in range(piters):
            nu = apply_op(gamma, u)
            r = nu.max()
            if r <= 0: return 0.0
            u = nu / r
        return r
    lo, hi = 0.4, 1.0
    for _ in range(18):
        mid = (lo + hi) / 2
        if rho(mid) >= 1.0: lo = mid
        else: hi = mid
    return lo

for K, T in [(3, 30), (4, 36), (5, 42)]:
    g = solve(K, T)
    out[f"K={K}"] = round(g, 4)
    print(f"K={K} T={T}: class-count gamma = {g:.4f}", flush=True)
print(json.dumps(out, indent=1))
