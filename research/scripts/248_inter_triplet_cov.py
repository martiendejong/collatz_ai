"""
248_inter_triplet_cov.py
========================
Verify the inter-triplet Cov structure revealed by Script 247.

For r=2 CODE-triplet group g (g equiv 2 mod 3):
  v2_at_sigma0 CODE-triplet = {v2[j*], v2[j*+Nl3], v2[j*+2Nl3]} where j* = 4g%Nl
  These are the three SLOTS of the v2 CODE-triplet at base h0 = j*%Nl3 (if j*<Nl3, slot=0).
  The cb CODE-triplet = {cb[j*], cb[j*+Nl3], cb[j*+2Nl3]}
  where cb[j*] = min(v2[h0], v2[h0+Nl3], v2[h0+2Nl3])           (min of triplet at base h0)
        cb[j*+Nl3] = min(v2[h0+Nl9], v2[h0+Nl9+Nl3], v2[h0+Nl9+2Nl3])  (min of triplet at h0+Nl9)
        cb[j*+2Nl3] = min(v2[h0+2Nl9], ...)                       (min of triplet at h0+2Nl9)
  (where Nl9 = Nl//9)

KEY REALIZATION (from Script 247):
  Slot 0 Cov: v2[j*] vs cb[j*]=min(v2[h0],..,v2[h0+2Nl3])
    v2[j*] IS one of the values in min(v2[h0],v2[h0+Nl3],v2[h0+2Nl3]) when j* = h0, h0+Nl3, or h0+2Nl3
    Actually for j* = 4g%Nl: j* is the v2 SLOT-0 value (if j*<Nl3) or SLOT-1 (if Nl3<=j*<2Nl3)
    or SLOT-2 (if 2Nl3<=j*<3Nl3). The slot of j* in the v2 CODE-triplet.
    cb[j*] uses the v2 CODE-triplet at base j*%Nl3 = j* - (j*//Nl3)*Nl3.
    
  CONCLUSION: Slot 0 Cov = INTRA-triplet Cov: v2 value vs min of OWN triplet.
    For i.i.d. exchangeable: Cov(X_slot0, min(X,Y,Z)) = (1/3)(E[min^2]-muE[min]) (same formula!).
    Observed POSITIVE (Script 246): strong positive within-triplet correlation dominates.

  Slot 1 Cov: v2[j*+Nl3] vs cb[j*+Nl3]=min of triplet at base h0+Nl9.
    v2[j*+Nl3] = SLOT-1 of v2 CODE-triplet at base j*%Nl3.
    cb[j*+Nl3] = min of DIFFERENT CODE-triplet at base j*%Nl3 + Nl9.
    This is an INTER-TRIPLET Cov between two DIFFERENT CODE-triplets.
    Observed NEGATIVE (-3.7e-3).

This script confirms the inter-triplet structure and checks if it holds universally.
"""
import numpy as np
from math import log2
import sys

ALPHA = log2(3.0)
N_ITER = 500

def run_kl(k, lam, n_iter=N_ITER):
    A  = lam ** -2.0
    B1 = lam ** (ALPHA - 2.0)
    B3 = lam ** (ALPHA - 1.0)
    N  = 3 ** (k - 1)
    Nl = N // 3
    i  = np.arange(N, dtype=np.int64)
    T4 = (4 * i + 2) % N
    s_arr, r_arr = np.divmod(i, 3)
    m0, m2 = (r_arr == 0), (r_arr == 2)
    R1 = (4 * s_arr) % Nl
    R3 = (2 * s_arr + 1) % Nl
    v = np.ones(N, dtype=np.float64)
    for _ in range(n_iter):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        w  = A * v[T4]
        w[m2] += B3 * cb[R3[m2]]
        w[m0] += B1 * cb[R1[m0]]
        v = w / w.max()
    return v, Nl

def analyze_inter_triplet(k, lam):
    v, Nl = run_kl(k, lam)
    v2 = v[2::3]  # length Nl
    cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
    
    Nl3 = Nl // 3  # CODE-triplet spacing in v2
    Nl9 = Nl // 9  # coarse base spacing (Nl3//3)
    
    # For r=2 groups g (g equiv 2 mod 3, g in [0, Nl3)):
    g_r2 = np.arange(Nl3, dtype=np.int64)[np.arange(Nl3)%3==2]
    j_star = (4 * g_r2) % Nl  # sigma0(g) for g equiv 2
    
    # The v2 CODE-triplet of v2_at_sigma0 at group g:
    # {v2[j*], v2[(j*+Nl3)%Nl], v2[(j*+2Nl3)%Nl]}
    v2_s0 = v2[j_star]
    v2_s1 = v2[(j_star + Nl3) % Nl]
    v2_s2 = v2[(j_star + 2*Nl3) % Nl]
    v2_mean = (v2_s0 + v2_s1 + v2_s2) / 3.0
    ld_v2_s0 = np.log2(v2_s0 / v2_mean)
    ld_v2_s1 = np.log2(v2_s1 / v2_mean)
    ld_v2_s2 = np.log2(v2_s2 / v2_mean)
    
    # The cb CODE-triplet at group g:
    # {cb[j*], cb[(j*+Nl3)%Nl], cb[(j*+2Nl3)%Nl]}
    cb_s0 = cb[j_star]
    cb_s1 = cb[(j_star + Nl3) % Nl]
    cb_s2 = cb[(j_star + 2*Nl3) % Nl]
    cb_mean = (cb_s0 + cb_s1 + cb_s2) / 3.0
    ld_cb_s0 = np.log2(cb_s0 / cb_mean)
    ld_cb_s1 = np.log2(cb_s1 / cb_mean)
    ld_cb_s2 = np.log2(cb_s2 / cb_mean)
    
    cov_s0 = float(np.mean(ld_v2_s0 * ld_cb_s0))
    cov_s1 = float(np.mean(ld_v2_s1 * ld_cb_s1))
    cov_s2 = float(np.mean(ld_v2_s2 * ld_cb_s2))
    cov_total = (cov_s0 + cov_s1 + cov_s2) / 3.0
    
    # INTRA-TRIPLET COV (slot 0): Cov(v2[j*], min(v2 at triplet of j*))
    # j* is in the v2 CODE-triplet at base h0 = j*%Nl3
    h0 = j_star % Nl3
    slot_of_jstar = j_star // Nl3
    # The v2 CODE-triplet of j* in v2_interleaved:
    triplet_v0 = v2[h0]
    triplet_v1 = v2[h0 + Nl3]
    triplet_v2 = v2[h0 + 2*Nl3]
    triplet_min = np.minimum(np.minimum(triplet_v0, triplet_v1), triplet_v2)
    # v2 value at j*: depending on slot, = v2[h0], v2[h0+Nl3], or v2[h0+2Nl3]
    v2_at_jstar = v2[j_star]  # = v2[h0 + slot*Nl3]
    mu_v2 = float(np.mean(v2_at_jstar))
    E_min = float(np.mean(triplet_min))
    # Intra-triplet Cov: Cov(v2[j*], min(v2-triplet-of-j*))
    cov_intra = float(np.mean((v2_at_jstar - mu_v2) * (triplet_min - E_min)))
    
    # INTER-TRIPLET COV (slot 1): Cov(v2[j*+Nl3], min(v2 at triplet of h0+Nl9))
    v2_slot1_val = v2[(j_star + Nl3) % Nl]
    h1 = (h0 + Nl9) % Nl3
    inter_triplet_min = np.minimum(np.minimum(v2[h1], v2[h1+Nl3]), v2[h1+2*Nl3])
    mu_v2s1 = float(np.mean(v2_slot1_val))
    E_inter_min = float(np.mean(inter_triplet_min))
    cov_inter = float(np.mean((v2_slot1_val - mu_v2s1) * (inter_triplet_min - E_inter_min)))
    
    return {
        'cov_s0': cov_s0, 'cov_s1': cov_s1, 'cov_s2': cov_s2, 'cov_total': cov_total,
        'cov_intra': cov_intra, 'cov_inter': cov_inter, 'Nl3': Nl3, 'Nl9': Nl9,
    }

print("248: Inter-triplet Cov structure for r=2 groups")
print("Confirm: slot 0 = intra-triplet Cov (pos), slot 1 = inter-triplet Cov (neg)")
print("="*80)

K = 8
LAM = 1.70
res = analyze_inter_triplet(K, LAM)
print(f"\nk={K}, lam={LAM}:")
print(f"  Slot 0 Cov: {res['cov_s0']:+.4e}")
print(f"  Slot 1 Cov: {res['cov_s1']:+.4e}")
print(f"  Slot 2 Cov: {res['cov_s2']:+.4e}")
print(f"  Total: {res['cov_total']:+.4e}")
print(f"  Intra-triplet Cov (raw): {res['cov_intra']:+.4e}")
print(f"  Inter-triplet Cov (slot1 vs shifted min, raw): {res['cov_inter']:+.4e}")

print()
print("Lambda scan at k=8:")
print(f"{'lam':>6}  {'cov_s0':>10}  {'cov_s1':>10}  {'cov_s2':>10}  {'total':>10}  {'cov_intra':>10}  {'cov_inter':>10}")
for lam in [1.30, 1.50, 1.70, 1.90, 2.00]:
    r = analyze_inter_triplet(8, lam)
    print(f"lam={lam:.2f}  {r['cov_s0']:>10.4e}  {r['cov_s1']:>10.4e}  {r['cov_s2']:>10.4e}  "
          f"{r['cov_total']:>10.4e}  {r['cov_intra']:>10.4e}  {r['cov_inter']:>10.4e}")
    sys.stdout.flush()

print()
print("Depth scan at lam=1.70:")
print(f"{'k':>4}  {'cov_s0':>10}  {'cov_s1':>10}  {'cov_s2':>10}  {'inter_neg?':>10}")
for k in range(5, 12):
    r = analyze_inter_triplet(k, 1.70)
    print(f"k={k:>2}  {r['cov_s0']:>10.4e}  {r['cov_s1']:>10.4e}  {r['cov_s2']:>10.4e}  "
          f"{'YES' if r['cov_inter']<0 else 'NO':>10}")
    sys.stdout.flush()

print()
print("done")
