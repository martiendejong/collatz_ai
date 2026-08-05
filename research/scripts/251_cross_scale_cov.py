"""
251_cross_scale_cov.py
======================
Corrected analysis of the slot-1 anti-correlation mechanism.

Script 247 shows: Cov(ld_v2_s1, ld_cb_s1) < 0 for r=2 groups, slots 1 and 2.
Script 250 was WRONG: it computed coarser-level intra-triplet Cov (positive by definition).

CORRECT UNDERSTANDING of slot1:
  v2_slot1 = v2[j*+Nl3]          where j* = (4g) % Nl, g≡2 mod 3
  cb_slot1 = cb[j*+Nl3]          = min(v2[m0], v2[m0+Nl3], v2[m0+2Nl3])
             where m0 = (j*+Nl3)//3 = j*//3 + Nl9

  KEY: j*+Nl3 = 3*m0 + 2  (exact, since j*≡2 mod 3 and Nl3≡0 mod 3)
  So v2[j*+Nl3] = v2[3*m0+2] — a "fine" v2 value at index 3*m0+2.
  And cb[j*+Nl3] = min(v2[m0], v2[m0+Nl3], v2[m0+2Nl3]) — the min at the "coarser" level m0.

CROSS-SCALE STRUCTURE:
  v2[3*m0+2] is at fine-level index 3*m0+2.
  cb uses v2 at coarser indices m0, m0+Nl3, m0+2Nl3.
  Note: 3*m0+2 is NOT in {m0, m0+Nl3, m0+2Nl3} in general.

THE QUESTION: why does large v2[3*m0+2] (relative to its sigma0 CODE-triplet)
correlate with small min(v2[m0], v2[m0+Nl3], v2[m0+2Nl3]) (relative to its meta-triplet)?

HYPOTHESIS: The K-L eigenvector has a SCALE-ALTERNATING structure:
  - At "fine" positions (3*m0+2): value ABOVE the sigma0 triplet mean
  => At "coarser" positions (m0): minimum BELOW the meta-triplet mean
  This would be a CROSS-SCALE ANTI-CORRELATION, consistent with the negative Cov.

Concretely: is v2[3*m0+2] and v2[m0] NEGATIVELY correlated?

Alternative test: is there a SIGN-FLIP between even and odd levels of the recursive structure?

Let p = 3*m0+2 (fine) and q = m0 (coarse). Note p = 3q+2, so q = p//3.
Hypothesis: v2[3q+2] and v2[q] are NEGATIVELY correlated within the eigenvector.

This would explain: ld of v2[3m0+2] up => v2[3m0+2] is large => v2[m0] is small
=> min(v2[m0], ...) is small => ld_cb < 0. QED.
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
    return v, Nl

def analyze_cross_scale(k, lam):
    v, Nl = run_kl(k, lam)
    v2 = v[2::3]   # v2[j] = v at (s=j, r=2), length Nl
    cb_arr = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
    Nl3 = Nl // 3
    Nl9 = Nl // 9

    # r=2 groups: g ≡ 2 mod 3, g ∈ [0, Nl3)
    g_r2 = np.arange(Nl3, dtype=np.int64)[np.arange(Nl3)%3 == 2]
    j_star = (4 * g_r2) % Nl   # sigma0 map

    # Verify j* ≡ 2 mod 3 always
    assert np.all(j_star % 3 == 2)

    # Slot 1: indices j*+Nl3 in v2 and cb
    js1 = (j_star + Nl3) % Nl
    assert np.all(js1 % 3 == 2)  # (j*+Nl3)%3 = (2+0)%3 = 2 (since Nl3%3=0)

    # v2 at slot1
    v2_s1 = v2[js1]   # v2 at fine index js1 = j*+Nl3 = 3*m0+2

    # cb at slot1
    cb_s1 = cb_arr[js1]   # = min(v2[m0], v2[m0+Nl3], v2[m0+2Nl3]) where m0 = js1//3

    m0_arr = js1 // 3   # m0 = (j*+Nl3)//3 = j*//3 + Nl9

    # Verify js1 = 3*m0 + 2 (holds when no Nl wrap-around in js1 = (j*+Nl3)%Nl)
    # For wrap-around cases: js1 = j*+Nl3-Nl, m0 = (j*+Nl3-Nl)//3 = j*//3+Nl9-Nl//3
    # The key property j*+Nl3 ≡ 2 mod 3 always (so js1 ≡ 2 mod 3), verified above

    # The three cb triplet positions: m0, m0+Nl3, m0+2Nl3
    v2_at_m0     = v2[m0_arr % Nl]
    v2_at_m0pNl3 = v2[(m0_arr + Nl3) % Nl]
    v2_at_m0p2Nl3= v2[(m0_arr + 2*Nl3) % Nl]

    # Verify cb_s1 = min of these
    cb_check = np.minimum(np.minimum(v2_at_m0, v2_at_m0pNl3), v2_at_m0p2Nl3)
    assert float(np.max(np.abs(cb_s1 - cb_check))) < 1e-10, "cb mismatch"

    # === HYPOTHESIS 1: Corr(v2[js1], v2[m0]) < 0 ===
    # v2[js1] = v2[3*m0+2] and v2[m0]: are these anti-correlated?
    corr_v2_js1_m0 = float(np.corrcoef(v2_s1, v2_at_m0)[0,1])
    corr_v2_js1_min = float(np.corrcoef(v2_s1, cb_s1)[0,1])

    print(f"\nk={k}, lam={lam}:")
    print(f"  Corr(v2[3m0+2], v2[m0])      = {corr_v2_js1_m0:+.4f}  (NEGATIVE => anti-corr hypothesis)")
    print(f"  Corr(v2[3m0+2], min_t(m0))   = {corr_v2_js1_min:+.4f}")
    print(f"  Cov(v2[3m0+2], v2[m0])       = {float(np.cov(v2_s1, v2_at_m0)[0,1]):+.4e}")
    print(f"  Cov(v2[3m0+2], min_t(m0))    = {float(np.cov(v2_s1, cb_s1)[0,1]):+.4e}")

    # === The log-deviation anti-correlation ===
    # sigma0-triplet CODE-variance ld:
    js0 = j_star
    js2 = (j_star + 2*Nl3) % Nl
    v2_s0 = v2[js0]; v2_s2 = v2[js2]
    v2_trip_mean = (v2_s0 + v2_s1 + v2_s2) / 3.0
    ld_v2_s1 = np.log2(v2_s1 / v2_trip_mean)

    cb_s0 = cb_arr[js0]; cb_s2 = cb_arr[js2]
    cb_trip_mean = (cb_s0 + cb_s1 + cb_s2) / 3.0
    ld_cb_s1 = np.log2(cb_s1 / cb_trip_mean)

    cov_ld = float(np.mean(ld_v2_s1 * ld_cb_s1))
    print(f"  Cov(ld_v2_s1, ld_cb_s1)      = {cov_ld:+.4e}  (NEGATIVE = anti-correlation, per Script 247)")

    # === Decompose: how does ld_v2_s1 relate to (v2[3m0+2] vs v2[m0])? ===
    # ld_v2_s1 > 0 means v2[3m0+2] > (v2[3m0_0+2] + v2[3m0+2] + v2[3m0_2+2])/3
    # where m0_0 = (j*+0)//3 = j*//3, m0_2 = (j*+2Nl3)//3 = j*//3 + 2Nl9

    m0_s0 = js0 // 3   # = j*//3
    m0_s2 = js2 // 3   # = j*//3 + 2Nl9

    v2_m0_s0 = v2[m0_s0 % Nl]
    v2_m0_s1 = v2_at_m0         # m0 for slot1
    v2_m0_s2 = v2[m0_s2 % Nl]

    # The three "m0" values form their OWN triplet in the coarser structure
    # (spacing Nl9 in m0-space):
    # m0_s0 = j*//3
    # m0_s1 = j*//3 + Nl9
    # m0_s2 = j*//3 + 2*Nl9
    # Note: m0_arr = js1//3; m0_s0 = js0//3 = j*//3. Due to modular wrap in js1,
    # m0_arr = m0_s0 + Nl9 only when j*+Nl3 < Nl. Otherwise m0_arr = m0_s0+Nl9-Nl//3.

    # These m0 values are spaced Nl9 apart: they form a coarse triplet
    coarse_trip_mean = (v2_m0_s0 + v2_m0_s1 + v2_m0_s2) / 3.0
    ld_v2_m0_s0 = np.log2(v2_m0_s0 / coarse_trip_mean)
    ld_v2_m0_s1 = np.log2(v2_m0_s1 / coarse_trip_mean)
    ld_v2_m0_s2 = np.log2(v2_m0_s2 / coarse_trip_mean)

    # Corr between fine ld (v2 at 3m0+2) and coarse ld (v2 at m0)
    corr_fine_coarse = float(np.corrcoef(ld_v2_s1, ld_v2_m0_s1)[0,1])
    print(f"\n  Corr(ld_v2[3m0+2 fine], ld_v2[m0 coarse]) = {corr_fine_coarse:+.4f}")
    print(f"  (if negative: cross-scale anti-correlation explains Cov<0)")

    # Now: how does ld_cb_s1 relate to ld_v2[m0]?
    # cb_s1 = min_t(m0) and cb_trip_mean = (min_t(m0_s0)+min_t(m0_s1)+min_t(m0_s2))/3
    # ld_cb_s1 = log(min_t(m0_s1) / coarse_min_mean)
    # min_t at coarse positions:
    min_t_m0_s0 = np.minimum(np.minimum(v2[m0_s0%Nl], v2[(m0_s0+Nl3)%Nl]), v2[(m0_s0+2*Nl3)%Nl])
    min_t_m0_s1 = cb_s1  # already computed
    min_t_m0_s2 = np.minimum(np.minimum(v2[m0_s2%Nl], v2[(m0_s2+Nl3)%Nl]), v2[(m0_s2+2*Nl3)%Nl])

    # Verify cb_trip_mean = (min_t_m0_s0 + min_t_m0_s1 + min_t_m0_s2)/3
    cb_trip_check = (min_t_m0_s0 + min_t_m0_s1 + min_t_m0_s2) / 3.0
    err = float(np.max(np.abs(cb_trip_mean - cb_trip_check)))
    print(f"  cb_trip_mean verification error: {err:.2e}  {'OK' if err<1e-8 else 'FAIL'}")

    # ld_cb_s1 in terms of min_t at coarse level:
    coarse_min_mean = (min_t_m0_s0 + min_t_m0_s1 + min_t_m0_s2) / 3.0
    ld_min_t_m0_s1 = np.log2(min_t_m0_s1 / coarse_min_mean)

    # This should equal ld_cb_s1
    assert float(np.max(np.abs(ld_min_t_m0_s1 - ld_cb_s1))) < 1e-8

    # Now: Cov(ld_v2[fine=3m0+2], ld_min_t[coarse=m0])
    # = Cov(ld_v2_s1, ld_cb_s1) = -3.7e-3 (from Script 247)

    # The FACTORED FORM: can we factor this as
    # Cov(ld_fine, ld_min_t_coarse) = Cov(ld_fine, ld_coarse_v2) * (some factor)?
    # where "ld_coarse_v2" is the ld of v2[m0] at the coarse level

    # Corr(ld_fine, ld_coarse_v2):
    print(f"\n  Corr(ld_v2_fine=3m0+2, ld_v2_coarse=m0) = {float(np.corrcoef(ld_v2_s1, ld_v2_m0_s1)[0,1]):+.4f}")
    print(f"  Corr(ld_v2_coarse=m0, ld_min_t_m0)     = {float(np.corrcoef(ld_v2_m0_s1, ld_min_t_m0_s1)[0,1]):+.4f}")
    print(f"  Corr(ld_v2_fine, ld_min_t_m0)           = {float(np.corrcoef(ld_v2_s1, ld_min_t_m0_s1)[0,1]):+.4f}")

    # Also: direct correlations v2[fine] vs v2[coarse] (raw, not log-dev)
    corr_raw = float(np.corrcoef(v2_s1, v2_at_m0)[0,1])
    print(f"\n  Corr(v2[3m0+2], v2[m0]) raw              = {corr_raw:+.4f}")

    # Summary
    print(f"\n  === SUMMARY ===")
    print(f"  v2[3m0+2] and v2[m0]: {'ANTI-CORR' if corr_v2_js1_m0 < 0 else 'POS-CORR'} ({corr_v2_js1_m0:+.4f})")
    print(f"  ld_v2_fine and ld_coarse_v2: {'ANTI-CORR' if corr_fine_coarse < 0 else 'POS-CORR'} ({corr_fine_coarse:+.4f})")
    return corr_v2_js1_m0, corr_fine_coarse, cov_ld

print("251: Cross-scale anti-correlation mechanism for slot-1 Cov<0")
print("Testing: does v2[3m0+2] anti-correlate with v2[m0]?")
print("="*70)

analyze_cross_scale(8, 1.70)

print("\n\n=== Lambda scan, k=8 ===")
print(f"{'lam':>6}  {'Corr_raw':>10}  {'Corr_ld':>10}  {'Cov_ld_s1':>12}")
for lam in [1.30, 1.50, 1.70, 1.90, 2.00]:
    v, Nl = run_kl(8, lam)
    v2 = v[2::3]
    cb_arr = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
    Nl3 = Nl//3; Nl9 = Nl//9
    g_r2 = np.arange(Nl3)[np.arange(Nl3)%3==2]
    j_star = (4*g_r2) % Nl
    js1 = (j_star + Nl3) % Nl
    m0_arr2 = js1 // 3

    v2_s1 = v2[js1]; cb_s1 = cb_arr[js1]
    v2_at_m0_2 = v2[m0_arr2 % Nl]
    corr_raw2 = float(np.corrcoef(v2_s1, v2_at_m0_2)[0,1])

    v2_s0 = v2[j_star]; v2_s2 = v2[(j_star+2*Nl3)%Nl]
    trip_mean = (v2_s0 + v2_s1 + v2_s2) / 3.0
    ld_v2_s1 = np.log2(v2_s1 / trip_mean)

    m0_s0 = j_star // 3; m0_s1 = m0_arr2; m0_s2 = (j_star+2*Nl3)//3
    v2_m0 = v2[m0_s1%Nl]; v2_m0_s0 = v2[m0_s0%Nl]; v2_m0_s2 = v2[m0_s2%Nl]
    coarse_mean = (v2_m0_s0 + v2_m0 + v2_m0_s2) / 3.0
    ld_v2_m0 = np.log2(v2_m0 / coarse_mean)
    corr_ld2 = float(np.corrcoef(ld_v2_s1, ld_v2_m0)[0,1])

    cb_s0 = cb_arr[j_star]; cb_s2 = cb_arr[(j_star+2*Nl3)%Nl]
    cb_mean = (cb_s0 + cb_s1 + cb_s2) / 3.0
    ld_cb_s1 = np.log2(cb_s1 / cb_mean)
    cov_ld2 = float(np.mean(ld_v2_s1 * ld_cb_s1))

    print(f"lam={lam:.2f}  {corr_raw2:>10.4f}  {corr_ld2:>10.4f}  {cov_ld2:>12.4e}")
    sys.stdout.flush()

print("\n\n=== Depth scan, lam=1.70 ===")
print(f"{'k':>4}  {'Corr_raw':>10}  {'Corr_ld':>10}  {'Cov_ld_s1':>12}")
LAM = 1.70
for k in range(5, 12):
    v, Nl = run_kl(k, LAM)
    v2 = v[2::3]
    cb_arr = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
    Nl3 = Nl//3; Nl9 = Nl//9
    g_r2 = np.arange(Nl3)[np.arange(Nl3)%3==2]
    j_star = (4*g_r2) % Nl
    js1 = (j_star+Nl3) % Nl
    m0_arr2 = js1 // 3
    v2_s1 = v2[js1]; cb_s1 = cb_arr[js1]
    v2_at_m0_2 = v2[m0_arr2 % Nl]
    corr_raw2 = float(np.corrcoef(v2_s1, v2_at_m0_2)[0,1])
    v2_s0 = v2[j_star]; v2_s2 = v2[(j_star+2*Nl3)%Nl]
    trip_mean = (v2_s0+v2_s1+v2_s2)/3.0
    ld_v2_s1 = np.log2(v2_s1/trip_mean)
    m0_s0=j_star//3; m0_s1=m0_arr2; m0_s2=(j_star+2*Nl3)//3
    v2_m0=v2[m0_s1%Nl]; v2_m0s0=v2[m0_s0%Nl]; v2_m0s2=v2[m0_s2%Nl]
    coarse_mean=(v2_m0s0+v2_m0+v2_m0s2)/3.0
    ld_v2_m0=np.log2(v2_m0/coarse_mean)
    corr_ld2=float(np.corrcoef(ld_v2_s1, ld_v2_m0)[0,1])
    cb_s0=cb_arr[j_star]; cb_s2=cb_arr[(j_star+2*Nl3)%Nl]
    cb_mean=(cb_s0+cb_s1+cb_s2)/3.0
    ld_cb_s1=np.log2(cb_s1/cb_mean)
    cov_ld2=float(np.mean(ld_v2_s1*ld_cb_s1))
    print(f"k={k:>2}  {corr_raw2:>10.4f}  {corr_ld2:>10.4f}  {cov_ld2:>12.4e}")
    sys.stdout.flush()

print("\ndone")
