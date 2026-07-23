"""
164_verify_KL_result.py
========================
VERIFICATION of the k>11 Krasikov-Lagarias exponents from script 163,
before any claim is made. Two independent checks:

 A. BRUTE-FORCE CROSS-CHECK of the compact index algebra: build the system
    directly from the paper's formulas (explicit residues m == 2 mod 3;
    branch by m mod 9; targets (4m-2)/3, (2m-1)/3 as integers reduced mod
    3^(k-1); cbar = min over the three lifts m', m'+3^(k-1), m'+2*3^(k-1)),
    with NO index tricks. Compare the Perron growth rho(lambda) against the
    compact implementation at several lambda for k = 3, 4, 5, 6.

 B. FEASIBILITY CERTIFICATE with explicit margin at k = 12 and k = 14:
    at lambda_cert (strictly below the bisection lambda*), converge the
    Perron vector v hard, rescale, and verify  F(v) >= (1 + delta) * v
    entrywise in float64, reporting the worst margin delta. A positive
    margin ~1e-3 dwarfs float rounding (~1e-15) and coefficient error,
    making the feasible LP solution solid at numerical-rigor level.
    (The remaining formal step for publication: interval arithmetic on the
    three lambda-power coefficients; noted, not done here.)
"""
import numpy as np
from math import log2

ALPHA = log2(3.0)

# ---------- compact implementation (as in 163) ----------
def make_maps(k):
    N = 3**(k-1)
    i = np.arange(N, dtype=np.int64)
    T4 = (4*i + 2) % N
    s, r = np.divmod(i, 3)
    m1 = (r == 0); m3 = (r == 2)
    Nl = N // 3
    R1 = np.where(m1, (4*s) % Nl, 0).astype(np.int64)
    R3 = np.where(m3, (2*s + 1) % Nl, 0).astype(np.int64)
    return N, Nl, T4, m1, m3, R1, R3

def step_compact(v, lmb, maps):
    N, Nl, T4, m1, m3, R1, R3 = maps
    cbar = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:3*Nl])
    w = v[T4] * lmb**(-2.0)
    w[m1] += cbar[R1[m1]] * lmb**(ALPHA-2.0)
    w[m3] += cbar[R3[m3]] * lmb**(ALPHA-1.0)
    return w

# ---------- brute-force from the paper ----------
def step_brute(vdict, lmb, k):
    """vdict: {m: value} for m in [3^k] (m == 2 mod 3). Direct formulas."""
    M3 = 3**k; M3l = 3**(k-1)
    def cbar(r):
        # r is a class mod 3^(k-1), r == 2 mod 3; lifts mod 3^k:
        return min(vdict[r % M3 if (r % M3) % 3 == 2 else r],
                   vdict[(r + M3l) % M3],
                   vdict[(r + 2*M3l) % M3])
    out = {}
    for m in vdict:
        w = vdict[(4*m) % M3] * lmb**(-2.0)
        if m % 9 == 2:
            t = ((4*m - 2)//3) % M3l
            w += cbar(t) * lmb**(ALPHA-2.0)
        elif m % 9 == 8:
            t = ((2*m - 1)//3) % M3l
            w += cbar(t) * lmb**(ALPHA-1.0)
        out[m] = w
    return out

print("A. BRUTE-FORCE CROSS-CHECK (compact vs direct-from-paper)")
rng = np.random.default_rng(7)
ok_all = True
for k in [3, 4, 5, 6]:
    maps = make_maps(k)
    N = maps[0]
    for lmb in [1.3, 1.6, 1.8]:
        v = rng.random(N) + 0.5
        w_c = step_compact(v, lmb, maps)
        vdict = {3*i+2: v[i] for i in range(N)}
        w_b = step_brute(vdict, lmb, k)
        w_b_arr = np.array([w_b[3*i+2] for i in range(N)])
        diff = float(np.max(np.abs(w_c - w_b_arr)))
        ok = diff < 1e-12
        ok_all &= ok
        print(f"  k={k} lambda={lmb}: max|compact-brute| = {diff:.2e}  {'OK' if ok else 'MISMATCH'}")
print(f"index algebra {'VERIFIED' if ok_all else 'BROKEN -- STOP'}")
print()

print("B. FEASIBILITY CERTIFICATES (margin check)")
for k, gamma_claim in [(12, 0.8531), (13, 0.8630), (14, 0.8724)]:
    maps = make_maps(k)
    N = maps[0]
    margin_lam = 2.0**(gamma_claim - 0.0008)     # certify slightly below claim
    v = np.ones(N)
    for it in range(3000):
        w = step_compact(v, margin_lam, maps)
        v = w / w.max()
    w = step_compact(v, margin_lam, maps)
    ratio = w / v
    worst = float(ratio.min())
    print(f"  k={k}: certify lambda={margin_lam:.6f} (gamma={log2(margin_lam):.4f}): "
          f"min F(v)/v = {worst:.6f}  -> {'FEASIBLE (certificate holds)' if worst >= 1.0 else 'not yet: increase iters or lower gamma'}")
print()
print("Interpretation: min F(v)/v >= 1 exhibits an explicit feasible solution of")
print("L_k^NT(lambda); by Krasikov-Lagarias Theorem 2.2 (published, peer-reviewed)")
print("this yields pi_1(x) > x^{log2 lambda} for all sufficiently large x.")
