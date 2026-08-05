"""
272_m2m_bias_structure.py
=========================
Two-track analysis of step (3b): CoV^2(v2) > CoV^2(v0) => m2m_v2 <= m2m_v0.

TRACK A: Mean bias (not max bias) from log-normal approximation.
  If E[m2m_v0 - m2m_LN_v0] > E[m2m_v2 - m2m_LN_v2] (i.e., bias_v0 > bias_v2
  in the mean), and gap_LN is positive, we have a *mean-level* proof.

TRACK B: Structural argument via B-coefficient ordering.
  v0 uses B1 = lam^(alpha-2), v2 uses B3 = lam^(alpha-1) = lam * B1.
  Claim: for a K-L column triplet (A*Z + B*W) where W = block_min(Z),
  larger B/A => smaller m2m. Test this directly by sweeping B/A.

TRACK C: Wide (k, lam) sweep for c2/c0 < R. Include k=5..14, lam=1.1..1.95.
"""
import numpy as np
from math import log2, sqrt, exp, pi
from scipy import stats, integrate
import sys

ALPHA = log2(3.0)

def run_kl(k, lam, n_iter=None):
    if n_iter is None:
        n_iter = 600 + 200*max(0, k-8)
    A  = lam**-2.0; B1 = lam**(ALPHA-2.0); B3 = lam**(ALPHA-1.0)
    N  = 3**(k-1); Nl = N//3
    i  = np.arange(N, dtype=np.int64)
    T4 = (4*i+2) % N
    s_arr, r_arr = np.divmod(i, 3)
    m0_mask, m2_mask = (r_arr==0), (r_arr==2)
    R1 = (4*s_arr) % Nl; R3 = (2*s_arr+1) % Nl
    v = np.ones(N, dtype=np.float64)
    for _ in range(n_iter):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        w  = A*v[T4]
        w[m2_mask] += B3*cb[R3[m2_mask]]
        w[m0_mask] += B1*cb[R1[m0_mask]]
        v = w/w.max()
    return v, Nl, A, B1, B3

def m2m_triplets(v0, Nl3):
    j3 = np.arange(Nl3, dtype=np.int64)
    col = np.stack([v0[j3], v0[j3+Nl3], v0[j3+2*Nl3]], axis=1)
    return float(np.mean(col.min(axis=1) / col.mean(axis=1)))

def get_c_R(v, Nl):
    Nl3 = Nl // 3
    v0 = v[0::3]; v1 = v[1::3]; v2 = v[2::3]
    j3 = np.arange(Nl3, dtype=np.int64)
    col0 = np.stack([v0[j3], v0[j3+Nl3], v0[j3+2*Nl3]], axis=1)
    col2 = np.stack([v2[j3], v2[j3+Nl3], v2[j3+2*Nl3]], axis=1)
    c0 = float(col0.min(axis=1).mean())
    c2 = float(col2.min(axis=1).mean())
    s = np.arange(Nl, dtype=np.int64)
    sigma1 = (4*s+2) % Nl
    lam_est = float(np.mean(v0)) / float(np.mean(v2)) * (float(np.mean(v1))/float(np.mean(v2)))
    A_val = v.max()**0  # A = lam^-2
    rho = A_val / float(np.mean(v1 / v0[sigma1]))  # Actually need A
    return c0, c2, float(np.mean(v0)), float(np.mean(v2))

def analyze_full(k, lam):
    A_val = lam**-2.0
    v, Nl, A, B1, B3 = run_kl(k, lam)
    Nl3 = Nl // 3
    v0 = v[0::3]; v1 = v[1::3]; v2 = v[2::3]
    j3 = np.arange(Nl3, dtype=np.int64)
    col0 = np.stack([v0[j3], v0[j3+Nl3], v0[j3+2*Nl3]], axis=1)
    col2 = np.stack([v2[j3], v2[j3+Nl3], v2[j3+2*Nl3]], axis=1)
    c0 = float(col0.min(axis=1).mean())
    c2 = float(col2.min(axis=1).mean())
    m0 = float(np.mean(v0)); m2 = float(np.mean(v2))
    s = np.arange(Nl, dtype=np.int64)
    sigma1 = (4*s+2) % Nl
    rho = A_val / float(np.mean(v1 / v0[sigma1]))
    t = A_val / rho
    R = (t**2 + lam) / (1 + t*lam)
    return c2/c0, R, c2/c0 <= R

# ======================================================================
# TRACK C: Wide sweep c2/c0 < R
# ======================================================================
print("="*68)
print("TRACK C: Wide sweep c2/c0 < R")
print(f"{'lam':>5}  k=5  k=6  k=7  k=8  k=10 k=12 k=14  all_ok")
lams = [1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 1.95]
for lam in lams:
    row = []
    ok_all = True
    for k in [5,6,7,8,10,12,14]:
        n_iter = 400 + 150*max(0, k-8)
        ratio, R, ok = analyze_full(k, lam)
        row.append("OK" if ok else "FAIL")
        if not ok:
            ok_all = False
    print(f"lam={lam:.2f}  " + "  ".join(f"{x:4s}" for x in row) + f"   {'ALL_OK' if ok_all else 'FAIL'}")
    sys.stdout.flush()

print()

# ======================================================================
# TRACK A: Mean signed bias
# ======================================================================
print("="*68)
print("TRACK A: Mean signed bias (LN approximation vs actual)")
print("bias = actual_m2m - LN_pred_m2m  (negative means LN overestimates)")
EZ1 = -0.84628  # E[min of 3 iid N(0,1)]

def track_a(k, lam):
    A_val = lam**-2.0
    v, Nl, A, B1, B3 = run_kl(k, lam)
    Nl3 = Nl // 3
    v0 = v[0::3]; v1 = v[1::3]; v2 = v[2::3]
    j3 = np.arange(Nl3, dtype=np.int64)
    col0 = np.stack([v0[j3], v0[j3+Nl3], v0[j3+2*Nl3]], axis=1)
    col2 = np.stack([v2[j3], v2[j3+Nl3], v2[j3+2*Nl3]], axis=1)
    def col_stats(col):
        lc = np.log(col + 1e-300)
        sig = lc.std(axis=1)
        mmr = col.min(axis=1) / col.mean(axis=1)
        ln_pred = np.exp(EZ1 * sig)
        bias = mmr - ln_pred
        return float(np.mean(mmr)), float(np.mean(ln_pred)), float(np.mean(bias)), float(np.std(bias))
    m2m0, ln0, bias0, std0 = col_stats(col0)
    m2m2, ln2, bias2, std2 = col_stats(col2)
    gap_actual = m2m0 - m2m2
    gap_ln = ln0 - ln2
    delta_bias = bias0 - bias2  # positive => bias_v0 > bias_v2 => actual gap > LN gap
    return gap_actual, gap_ln, delta_bias, bias0, bias2

print(f"{'lam':>5} {'k':>3}  gap_actual  gap_LN  delta_bias  bias0   bias2  proof_ok")
for lam in [1.30, 1.50, 1.70, 1.90]:
    for k in [8, 10, 12]:
        g_act, g_ln, d_bias, b0, b2 = track_a(k, lam)
        # proof: if gap_LN + delta_bias > 0 (which = gap_actual)
        ok = g_act > 0
        print(f"lam={lam:.2f} k={k:>2}  {g_act:>9.5f}  {g_ln:>7.5f}  {d_bias:>9.5f}  {b0:>7.4f}  {b2:>7.4f}  {ok}")
    sys.stdout.flush()

print()

# ======================================================================
# TRACK B: B-ratio structural argument
# Model: column triplet X_r = (A*Z_r + B*W_r) / rho
# where Z_r ~ some distribution and W_r = min of column of Z.
# Test: does larger B/A => smaller m2m?
# ======================================================================
print("="*68)
print("TRACK B: B/A ratio structural argument")
print("Synthetic: X = (A*Z + B*W) / rho where W = min(col of Z)")
print("Does larger B/A => smaller m2m?")

# Generate synthetic data: Z ~ uniform or log-normal within column triplets
rng = np.random.default_rng(12345)
N_cols = 10000

def synth_m2m(BA_ratio, sigma=0.2, seed=42):
    """Compute m2m for synthetic triplets X_r = A*Z_r + B*W_r
    where Z_r are correlated within column (block min W = min Z).
    Z columns are iid log-normal (0, sigma).
    B/A = BA_ratio (A=1 wlog, B=BA_ratio).
    """
    rng2 = np.random.default_rng(seed)
    # Generate 3 columns of Z (shape N_cols x 3)
    Z = np.exp(rng2.normal(0, sigma, size=(N_cols, 3)))
    W = Z.min(axis=1, keepdims=True)  # block min per column
    B = BA_ratio
    X = (1.0 * Z + B * W)
    X_mean = X.mean(axis=1)
    X_min = X.min(axis=1)
    return float(np.mean(X_min / X_mean))

# B1/A and B3/A for several lambda values
print(f"{'lam':>5}  B1/A   B3/A  m2m(B1) m2m(B3)  B3>B1 => m2m3<m2m1?")
for lam in [1.1, 1.3, 1.5, 1.7, 1.9]:
    A_val = lam**-2
    B1_val = lam**(ALPHA-2)
    B3_val = lam**(ALPHA-1)
    BA1 = B1_val / A_val
    BA3 = B3_val / A_val
    m0 = synth_m2m(BA1, sigma=0.10)
    m2 = synth_m2m(BA3, sigma=0.10)
    ok = m2 < m0
    print(f"lam={lam:.1f}  {BA1:.4f} {BA3:.4f}  {m0:.6f} {m2:.6f}  {ok}")
sys.stdout.flush()

print("\nB/A sweep (fixed sigma=0.10):")
print(f"{'B/A':>6}  m2m")
for ba in [0.0, 0.1, 0.3, 0.5, 0.8, 1.0, 1.5, 2.0, 3.0, 5.0]:
    m = synth_m2m(ba, sigma=0.10)
    print(f"{ba:>6.2f}  {m:.6f}")

print("\ndone")
