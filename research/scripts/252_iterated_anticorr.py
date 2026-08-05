"""
252_iterated_anticorr.py
========================
Test: is the parent-child anti-correlation Corr(v2[3s+2], v2[s]) < 0 a
ONE-LEVEL or MULTI-LEVEL phenomenon?

If Corr(v2[3s+2], v2[s]) < 0 (one level, proved empirically in Script 251),
does Corr(v2[9s+8], v2[3s+2]) < 0 as well? (iterating 3s+2 -> 3(3s+2)+2 = 9s+8)
And: Corr(v2[9s+8], v2[s]) > 0? (two anti-correlations compose to give positive)

Also: decompose the K-L equation to identify WHICH TERM drives the anti-correlation:
  ρ v2[s]     = A v[T4(3s+2)]   + B3 cb[(2s+1)%Nl]  ... (P1)
  ρ v2[3s+2]  = A v[T4(9s+8)]   + B3 cb[(6s+5)%Nl]  ... (P2)

  Anti-corr source: T4 term, cb term, or both?
"""
import numpy as np
from math import log2
import sys

ALPHA = log2(3.0)
N_ITER = 600

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
    return v, Nl, T4, R1, R3, A, B1, B3, N

def analyze_iterated(k, lam):
    v, Nl, T4, R1, R3, A, B1, B3, N = run_kl(k, lam)
    v2 = v[2::3]  # v2[s] = v at (s, r=2), length Nl
    cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
    Nl3 = Nl // 3
    rho = float(np.max(v))  # max normalization, so rho = 1 by construction. Actual rho from eq.

    # === ITERATED ANTI-CORRELATION ===
    # Level 0 -> 1: Corr(v2[3s+2], v2[s]) for s in [0, Nl3) (to ensure 3s+2 < Nl)
    s0 = np.arange(Nl//9, dtype=np.int64)  # small range to ensure valid indices
    # Actually, just use all s where 3s+2 < Nl and 9s+8 < Nl:
    max_s = (Nl - 9) // 9  # 9s+8 < Nl => s < (Nl-8)/9
    s0 = np.arange(max_s, dtype=np.int64)

    v2_s     = v2[s0]
    v2_3sp2  = v2[3*s0 + 2]        # level 1: 3s+2
    v2_9sp8  = v2[9*s0 + 8]        # level 2: 9s+8 = 3(3s+2)+2

    corr_01 = float(np.corrcoef(v2_3sp2, v2_s)[0,1])       # level 1 vs level 0
    corr_12 = float(np.corrcoef(v2_9sp8, v2_3sp2)[0,1])    # level 2 vs level 1
    corr_02 = float(np.corrcoef(v2_9sp8, v2_s)[0,1])       # level 2 vs level 0

    print(f"\nk={k}, lam={lam}, max_s={max_s}:")
    print(f"  Corr(v2[3s+2], v2[s])       = {corr_01:+.4f}  (level 1 vs 0)")
    print(f"  Corr(v2[9s+8], v2[3s+2])    = {corr_12:+.4f}  (level 2 vs 1)")
    print(f"  Corr(v2[9s+8], v2[s])       = {corr_02:+.4f}  (level 2 vs 0; POSITIVE if alternating)")

    # Check alternating sign: neg * neg = pos?
    print(f"  Product corr_01 * corr_12    = {corr_01*corr_12:+.4f}  (expected same sign as corr_02)")

    # === DECOMPOSE WHICH TERM DRIVES ANTI-CORR ===
    # K-L equation for v2[s] (s in s0):
    # i(s) = 3s+2 (the index into v). T4(i) = (4i+2)%N = (12s+10)%N.
    i_s = 3*s0 + 2
    T4_s = (4*i_s + 2) % N   # T4(3s+2) = (12s+10)%N

    # K-L equation for v2[3s+2]:
    i_3sp2 = 3*(3*s0+2) + 2  # = 9s+8
    T4_3sp2 = (4*i_3sp2 + 2) % N  # = (36s+34)%N

    # R3 for s: R3(s) = (2s+1)%Nl
    R3_s = (2*s0 + 1) % Nl
    # R3 for 3s+2: R3(3s+2) = (2*(3s+2)+1)%Nl = (6s+5)%Nl
    R3_3sp2 = (6*s0 + 5) % Nl

    # T4 pullback values
    v_T4_s     = v[T4_s]    # A*v[T4(3s+2)]
    v_T4_3sp2  = v[T4_3sp2] # A*v[T4(9s+8)]

    # cb values
    cb_s     = cb[R3_s]     # B3*cb[(2s+1)]
    cb_3sp2  = cb[R3_3sp2]  # B3*cb[(6s+5)]

    # Reconstruction (verify K-L equation)
    v2_s_recon = A * v_T4_s + B3 * cb_s
    v2_3sp2_recon = A * v_T4_3sp2 + B3 * cb_3sp2
    # Note: these give ρ*v2[s], not v2[s]. But ρ ≈ 1 (normalization). Close enough.

    # Correlations of COMPONENTS:
    corr_T4_T4 = float(np.corrcoef(v_T4_s, v_T4_3sp2)[0,1])
    corr_cb_cb = float(np.corrcoef(cb_s, cb_3sp2)[0,1])
    corr_T4_v2s = float(np.corrcoef(v_T4_s, v2_s)[0,1])
    corr_cb_v2s = float(np.corrcoef(cb_s, v2_s)[0,1])
    corr_T4_v23sp2 = float(np.corrcoef(v_T4_3sp2, v2_3sp2)[0,1])
    corr_cb_v23sp2 = float(np.corrcoef(cb_3sp2, v2_3sp2)[0,1])

    print(f"\n  === Decomposition of anti-correlation ===")
    print(f"  Corr(v[T4(3s+2)],  v[T4(9s+8)]) = {corr_T4_T4:+.4f}  (T4 terms corr)")
    print(f"  Corr(cb[(2s+1)],   cb[(6s+5)])   = {corr_cb_cb:+.4f}  (cb terms corr)")
    print(f"")
    print(f"  Corr(v2[s], v[T4(3s+2)])          = {corr_T4_v2s:+.4f}  (T4 input vs v2[s])")
    print(f"  Corr(v2[s], cb[(2s+1)])            = {corr_cb_v2s:+.4f}  (cb input vs v2[s])")
    print(f"")
    print(f"  Corr(v2[3s+2], v[T4(9s+8)])       = {corr_T4_v23sp2:+.4f}  (T4 input vs v2[3s+2])")
    print(f"  Corr(v2[3s+2], cb[(6s+5)])         = {corr_cb_v23sp2:+.4f}  (cb input vs v2[3s+2])")

    # Cross-correlations: does T4 of level 1 anti-correlate with cb of level 0?
    corr_T4_3sp2_cb_s = float(np.corrcoef(v_T4_3sp2, cb_s)[0,1])
    corr_cb_3sp2_T4_s = float(np.corrcoef(cb_3sp2, v_T4_s)[0,1])
    corr_T4_3sp2_v2_s = float(np.corrcoef(v_T4_3sp2, v2_s)[0,1])
    corr_cb_3sp2_v2_s = float(np.corrcoef(cb_3sp2, v2_s)[0,1])

    print(f"\n  Cross-term anti-corrs (WHICH DRIVES Corr(v2[3s+2], v2[s]) < 0):")
    print(f"  Corr(v2[s], v[T4(9s+8)])          = {corr_T4_3sp2_v2_s:+.4f}  (T4 of child vs parent)")
    print(f"  Corr(v2[s], cb[(6s+5)])            = {corr_cb_3sp2_v2_s:+.4f}  (cb of child vs parent)")
    print(f"  Corr(v[T4(9s+8)], cb[(2s+1)])     = {corr_T4_3sp2_cb_s:+.4f}")
    print(f"  Corr(cb[(6s+5)], v[T4(3s+2)])     = {corr_cb_3sp2_T4_s:+.4f}")

    # cb[(6s+5)] = min(v2[2s+1], v2[2s+1+Nl3], v2[2s+1+2Nl3]) — uses v2 at ~2s+1
    # cb[(2s+1)] uses v2 at ~(2s+1)/3
    # Direct: Corr(v2[2s+1], v2[s]) — is NEARBY v2 anti-correlated with parent?
    v2_2sp1 = v2[(2*s0+1) % Nl]
    v2_2sp1_div3 = v2[((2*s0+1)//3) % Nl]
    corr_2sp1_s = float(np.corrcoef(v2_2sp1, v2_s)[0,1])
    corr_2sp1d3_s = float(np.corrcoef(v2_2sp1_div3, v2_s)[0,1])
    print(f"\n  Corr(v2[2s+1], v2[s])             = {corr_2sp1_s:+.4f}  (cb-of-child neighbor vs parent)")
    print(f"  Corr(v2[(2s+1)/3], v2[s])         = {corr_2sp1d3_s:+.4f}  (cb-of-parent neighbor vs parent)")

    return corr_01, corr_12, corr_02

print("252: Iterated parent-child anti-correlation and decomposition")
print("="*70)

analyze_iterated(8, 1.70)

print("\n\n=== Lambda scan, k=8 ===")
print(f"{'lam':>6}  {'corr_01':>10}  {'corr_12':>10}  {'corr_02':>10}  {'prod01x12':>12}")
for lam in [1.30, 1.50, 1.70, 1.90, 2.00]:
    v, Nl, *_ = run_kl(8, lam)
    v2 = v[2::3]
    max_s = (Nl - 9) // 9
    s0 = np.arange(max_s, dtype=np.int64)
    c01 = float(np.corrcoef(v2[3*s0+2], v2[s0])[0,1])
    c12 = float(np.corrcoef(v2[9*s0+8], v2[3*s0+2])[0,1])
    c02 = float(np.corrcoef(v2[9*s0+8], v2[s0])[0,1])
    print(f"lam={lam:.2f}  {c01:>10.4f}  {c12:>10.4f}  {c02:>10.4f}  {c01*c12:>12.4f}")
    sys.stdout.flush()

print("\n\n=== Depth scan, lam=1.70 ===")
print(f"{'k':>4}  {'corr_01':>10}  {'corr_12':>10}  {'corr_02':>10}")
for k in range(6, 13):
    v, Nl, *_ = run_kl(k, 1.70)
    v2 = v[2::3]
    max_s = (Nl - 9) // 9
    s0 = np.arange(max(max_s, 1), dtype=np.int64)
    c01 = float(np.corrcoef(v2[3*s0+2], v2[s0])[0,1])
    c12 = float(np.corrcoef(v2[9*s0+8], v2[3*s0+2])[0,1])
    c02 = float(np.corrcoef(v2[9*s0+8], v2[s0])[0,1])
    print(f"k={k:>2}  {c01:>10.4f}  {c12:>10.4f}  {c02:>10.4f}")
    sys.stdout.flush()

print("\ndone")
