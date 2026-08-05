"""
250_ordering_analysis.py
========================
Investigate WHY the inter-triplet Cov is negative (slot 1,2 in Script 247).

Core question: when v2[j*+Nl3] is ABOVE its CODE-triplet mean (positive ld),
why does cb at position j*+Nl3 tend to be BELOW its CODE-triplet mean (negative ld)?

Key structural analysis:
- j* = (4g) % Nl for g ≡ 2 mod 3
- j* % 3 = 2 always (since 4g ≡ g mod 3, g ≡ 2 mod 3)
- cb[j] where j = 3s+r groups v-values by r-type over s-CODE-triplet
- v2[j*+Nl3] = v at (s=s_prime+Nl9, r=2) where s_prime = j*//3
- cb[j*+Nl3] = min of v at (s_prime+Nl9, r=2), (s_prime+Nl9+Nl3, r=2), (s_prime+Nl9+2Nl3, r=2)
  => v2[j*+Nl3] IS in the cb[j*+Nl3] min triplet (slot 0 of that triplet)!
  => This is Cov(X, min(X,Y,Z)) -- should be POSITIVE (Script 246)!
  => But log-deviation Cov is NEGATIVE. WHY?

Resolution: The log-deviation is taken RELATIVE TO THE sigma0-CODE-TRIPLET MEAN.
  sigma0-triplet of slot1: {j*+Nl3, j*+2Nl3, j*+3Nl3%Nl}
  cb-triplet of slot1: {cb[j*+Nl3], cb[j*+2Nl3], cb[j*+3Nl3%Nl]}
  These cb values are mins of DIFFERENT s-CODE-triplets (NOT the sigma0 triplet)!

  The anti-correlation arises from the interaction between:
  - ld_v2 = log(v2[j*+Nl3] / mean_sigma0_triplet) = slot0 within v2-meta-triplet
  - ld_cb = log(cb[j*+Nl3] / mean_cb_meta-triplet) = slot0 within cb-meta-triplet

  These DIFFER because the cb values at slots 0,1,2 of the sigma0-triplet
  are mins of DIFFERENT s-CODE-triplets, creating cross-triplet structure.
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

def analyze(k, lam):
    v, Nl = run_kl(k, lam)
    v2 = v[2::3]  # v2[s] = v at (s, r=2), length Nl
    cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
    Nl3 = Nl // 3
    Nl9 = Nl // 9

    print(f"\nk={k}, lam={lam}, Nl={Nl}, Nl3={Nl3}, Nl9={Nl9}")

    # === ORDERING within v2 CODE-triplets (s-triplets) ===
    h0 = np.arange(Nl3, dtype=np.int64)
    a0 = v2[h0]
    a1 = v2[h0 + Nl3]
    a2 = v2[h0 + 2*Nl3]
    print("\n=== v2 s-CODE-triplet ordering ===")
    print(f"  P(v2[s0] > v2[s1]) = {float(np.mean(a0 > a1)):.3f}")
    print(f"  P(v2[s1] > v2[s2]) = {float(np.mean(a1 > a2)):.3f}")
    print(f"  P(v2[s0] > v2[s2]) = {float(np.mean(a0 > a2)):.3f}")
    print(f"  slot0=max: {float(np.mean((a0>=a1)&(a0>=a2))):.3f}")
    print(f"  slot1=max: {float(np.mean((a1>=a0)&(a1>=a2))):.3f}")
    print(f"  slot2=max: {float(np.mean((a2>=a0)&(a2>=a1))):.3f}")

    # === r=2 group analysis (core of Cov mechanism) ===
    g_r2 = np.arange(Nl3, dtype=np.int64)[np.arange(Nl3)%3 == 2]
    j_star = (4 * g_r2) % Nl

    # Verify j* % 3 == 2 always
    assert np.all(j_star % 3 == 2), "j* should always be ≡ 2 mod 3"

    # s_prime = j*//3 (the s-position of j* in v2/cb indexing)
    s_prime = j_star // 3  # s' for j*

    # Slot 1: j*+Nl3
    js1 = (j_star + Nl3) % Nl
    s_js1 = js1 // 3  # s' for j*+Nl3 = s_prime + Nl9

    # Verify s_js1 = s_prime + Nl9 (mod Nl//3)
    Nl_div3 = Nl // 3
    assert np.all(s_js1 % Nl_div3 == (s_prime + Nl9) % Nl_div3), "s_js1 should be s_prime + Nl9"

    # v2[j*+Nl3] = v2[s_js1] = v at (s_prime+Nl9, r=2)
    # cb[j*+Nl3]: j*+Nl3 = 3s_prime*1+r2+Nl3... let's compute directly
    # cb[j*+Nl3] uses the s-CODE-triplet: {s_js1, s_js1+Nl3, s_js1+2Nl3} at r=2
    # min_t_at_s_js1 = min(v2[s_js1], v2[s_js1+Nl3], v2[s_js1+2*Nl3])

    # But s_js1 can be up to Nl//3-1 + Nl9. We need modular indexing.
    # Actually v2 has length Nl = N//3. s goes from 0 to Nl-1.
    # The s-triplet of j*+Nl3 is: s_js1, s_js1+Nl3, s_js1+2Nl3 (all mod Nl, but for cb we use cb[j])
    # cb[j*+Nl3] is already computed directly as cb[(j_star+Nl3)%Nl]

    cb_s1_direct = cb[(j_star + Nl3) % Nl]
    # This should equal min(v2[s_js1], v2[(s_js1+Nl3)%Nl], v2[(s_js1+2*Nl3)%Nl]) -- AT THE SAME r-TYPE
    # Actually: cb[j] = min(v[j], v[j+Nl], v[j+2Nl]) where v is the full 3Nl-vector
    # v[j*+Nl3] = v at i=(j*+Nl3). r-type = (j*+Nl3)%3 = (j*%3 + Nl3%3)%3 = (2 + 0)%3 = 2
    # (since Nl3 = Nl//3 = N//9 = 3^(k-3), divisible by 3)
    # s-type = (j*+Nl3)//3 = j*//3 + Nl3//3 = s_prime + Nl9
    # So cb[j*+Nl3] = min(v at (s_prime+Nl9, r=2), v at (s_prime+Nl9+Nl//3, r=2), v at (s_prime+Nl9+2Nl//3, r=2))
    # = min(v2[s_prime+Nl9], v2[s_prime+Nl9+Nl3], v2[s_prime+Nl9+2Nl3])
    # = min_t(s_prime+Nl9)

    # Verify this:
    min_t_check = np.minimum(np.minimum(
        v2[s_js1 % (Nl)],
        v2[(s_js1 + Nl3) % (Nl)],
    ), v2[(s_js1 + 2*Nl3) % (Nl)])
    # (taking mod Nl for safety)
    print(f"\n=== Verification: cb[j*+Nl3] == min_t(s_js1) ===")
    max_err = float(np.max(np.abs(cb_s1_direct - min_t_check)))
    print(f"  max |cb[j*+Nl3] - min_t(s_js1)| = {max_err:.2e}  {'OK' if max_err < 1e-10 else 'FAIL'}")

    # === THE KEY RELATIONSHIP ===
    # ld_v2_s1 = log(v2[s_prime+Nl9] / mean(v2[s_prime+Nl9], v2[s_prime+2Nl9], v2[s_prime+3Nl9]))
    # ld_cb_s1 = log(min_t(s_prime+Nl9) / mean(min_t(s_prime+Nl9), min_t(s_prime+2Nl9), min_t(s_prime+3Nl9)))
    # Here: s_prime+Nl9 = s_js1, s_prime+2Nl9 = s_js1+Nl9, s_prime+3Nl9 = s_js1+2Nl9

    # sigma0-triplet slots:
    # slot0: j*, slot1: j*+Nl3, slot2: j*+2Nl3
    # s_positions: s_prime, s_prime+Nl9, s_prime+2Nl9

    s0 = s_prime
    s1 = (s_prime + Nl9) % Nl
    s2 = (s_prime + 2*Nl9) % Nl

    # v2 at sigma0 slots
    v2_ss0 = v2[s0]
    v2_ss1 = v2[s1]
    v2_ss2 = v2[s2]
    v2_meta_mean = (v2_ss0 + v2_ss1 + v2_ss2) / 3.0
    ld_v2_s0 = np.log2(v2_ss0 / v2_meta_mean)
    ld_v2_s1_comp = np.log2(v2_ss1 / v2_meta_mean)
    ld_v2_s2_comp = np.log2(v2_ss2 / v2_meta_mean)

    # min_t at sigma0 slots (min of v2-triplet at each s-position)
    mt_s0 = np.minimum(np.minimum(v2[s0], v2[(s0+Nl3)%Nl]), v2[(s0+2*Nl3)%Nl])
    mt_s1 = np.minimum(np.minimum(v2[s1], v2[(s1+Nl3)%Nl]), v2[(s1+2*Nl3)%Nl])
    mt_s2 = np.minimum(np.minimum(v2[s2], v2[(s2+Nl3)%Nl]), v2[(s2+2*Nl3)%Nl])
    mt_meta_mean = (mt_s0 + mt_s1 + mt_s2) / 3.0
    ld_mt_s0 = np.log2(mt_s0 / mt_meta_mean)
    ld_mt_s1 = np.log2(mt_s1 / mt_meta_mean)
    ld_mt_s2 = np.log2(mt_s2 / mt_meta_mean)

    print(f"\n=== Core covariances ===")
    print(f"  Cov(ld_v2_s1, ld_mt_s1) = {float(np.mean(ld_v2_s1_comp * ld_mt_s1)):.4e}")
    print(f"  Cov(ld_v2_s0, ld_mt_s0) = {float(np.mean(ld_v2_s0 * ld_mt_s0)):.4e}  (slot0)")
    print(f"  Cov(ld_v2_s2, ld_mt_s2) = {float(np.mean(ld_v2_s2_comp * ld_mt_s2)):.4e}  (slot2)")

    # Now: decompose why ld_v2_s1 and ld_mt_s1 anti-correlate
    # ld_v2_s1 high means v2[s1] is large relative to v2 at (s0, s1, s2)
    # ld_mt_s1 high means min_t(s1) is large relative to min_t at (s0, s1, s2)
    # min_t(s) = min(v2[s], v2[s+Nl3], v2[s+2Nl3])

    # If v2[s1] is large (ld_v2_s1 > 0), what about v2[s1+Nl3] and v2[s1+2Nl3]?
    # These are at positions s1+Nl3, s1+2Nl3 -- NOT in the meta-triplet (s0, s1, s2).
    # The meta-triplet is spaced by Nl9, while the v2-triplet is spaced by Nl3 = 3*Nl9.

    # Key: ORDERING within the v2-triplet at s1:
    # v2[s1] vs v2[s1+Nl3] vs v2[s1+2Nl3]
    # How does ld_v2_s1 correlate with these?

    v2_s1_p1 = v2[(s1+Nl3)%Nl]   # v2 at s1+Nl3 (triplet neighbor)
    v2_s1_p2 = v2[(s1+2*Nl3)%Nl]  # v2 at s1+2Nl3 (triplet neighbor)

    # ld_v2_s1 vs v2 at triplet neighbors
    print(f"\n=== Cross-level structure at s1 ===")
    print(f"  Corr(ld_v2_s1, v2[s1+Nl3]) = {float(np.corrcoef(ld_v2_s1_comp, v2_s1_p1)[0,1]):.4f}")
    print(f"  Corr(ld_v2_s1, v2[s1+2Nl3]) = {float(np.corrcoef(ld_v2_s1_comp, v2_s1_p2)[0,1]):.4f}")

    # The min_t(s1) = min(v2[s1], v2[s1+Nl3], v2[s1+2Nl3])
    # If ld_v2_s1 > 0 (v2[s1] is above meta-mean), is min_t(s1) below meta-mean?
    # Only if v2[s1+Nl3] or v2[s1+2Nl3] are small enough to make min_t(s1) small.
    # But if ld_v2_s1 > 0, v2[s1] itself is large, so it's NOT the minimum.
    # min_t(s1) = min(v2[s1], v2[s1+Nl3], v2[s1+2Nl3]) <= v2[s1+Nl3] and <= v2[s1+2Nl3]

    # If v2[s1] is large (v2[s1] > meta-mean), min_t(s1) = min(v2[s1+Nl3], v2[s1+2Nl3]) ≈ min of the OTHER TWO!
    # So: ld_v2_s1 high => v2[s1] is NOT the min => min_t(s1) = min(v2[s1+Nl3], v2[s1+2Nl3])
    # Now, min_t(s1) RELATIVE TO meta-mean(min_t(s0), min_t(s1), min_t(s2)):
    # This depends on the cross-level ordering of v2 values.

    # Simpler test: Cov(v2[s1], min_t(s1)) vs Cov(ld_v2_s1, ld_mt_s1)
    cov_raw = float(np.cov(v2_ss1, mt_s1)[0,1])
    print(f"\n  Cov(v2[s1], min_t(s1)) = {cov_raw:.4e}  {'POS' if cov_raw>0 else 'NEG'}")
    print(f"  Cov(ld_v2_s1, ld_mt_s1) = {float(np.mean(ld_v2_s1_comp * ld_mt_s1)):.4e}")

    # The sign reversal from raw to ld is the KEY PHENOMENON
    # Let me check what sign Cov(X, min(X,Y,Z)) has vs Cov(ld_X, ld_min)
    # For the SAME triplet: Cov(v2[s1], min_t(s1)) should be POSITIVE (X in min triplet)
    # But ld Cov is NEGATIVE

    # Conditional test: given ld_v2_s1 > 0, what is E[ld_mt_s1]?
    pos_mask = ld_v2_s1_comp > 0
    neg_mask = ~pos_mask
    print(f"\n=== Conditional expectations ===")
    print(f"  E[ld_mt_s1 | ld_v2_s1 > 0] = {float(np.mean(ld_mt_s1[pos_mask])):.4e}")
    print(f"  E[ld_mt_s1 | ld_v2_s1 < 0] = {float(np.mean(ld_mt_s1[neg_mask])):.4e}")
    print(f"  (positive ld_v2 => NEGATIVE ld_mt? confirms anti-correlation)")

    # ORDERING: when v2[s1] is the largest in meta-triplet, what is min_t(s1) relative to meta-mean?
    s1_is_largest = (v2_ss1 >= v2_ss0) & (v2_ss1 >= v2_ss2)
    s1_is_smallest = (v2_ss1 <= v2_ss0) & (v2_ss1 <= v2_ss2)
    print(f"\n=== min_t ordering when v2[s1] is extreme ===")
    if s1_is_largest.sum() > 0:
        print(f"  E[ld_mt_s1 | v2[s1]=max in meta] = {float(np.mean(ld_mt_s1[s1_is_largest])):.4e}")
    if s1_is_smallest.sum() > 0:
        print(f"  E[ld_mt_s1 | v2[s1]=min in meta] = {float(np.mean(ld_mt_s1[s1_is_smallest])):.4e}")

# Run main analysis
analyze(8, 1.70)

print("\n\n=== Lambda scan: sign of Cov(v2[s1], min_t(s1)) ===")
for lam in [1.30, 1.50, 1.70, 1.90, 2.00]:
    v, Nl = run_kl(8, lam)
    v2 = v[2::3]
    Nl3, Nl9 = Nl//3, Nl//9
    g_r2 = np.arange(Nl3)[np.arange(Nl3)%3==2]
    j_star = (4*g_r2) % Nl
    s_prime = j_star // 3
    s1 = (s_prime + Nl9) % Nl
    v2_s1 = v2[s1]
    mt_s1 = np.minimum(np.minimum(v2[s1], v2[(s1+Nl3)%Nl]), v2[(s1+2*Nl3)%Nl])
    cov_raw = float(np.cov(v2_s1, mt_s1)[0,1])
    # meta triplet
    s0 = s_prime
    s2 = (s_prime + 2*Nl9) % Nl
    mt_s0 = np.minimum(np.minimum(v2[s0], v2[(s0+Nl3)%Nl]), v2[(s0+2*Nl3)%Nl])
    mt_s2 = np.minimum(np.minimum(v2[s2], v2[(s2+Nl3)%Nl]), v2[(s2+2*Nl3)%Nl])
    mt_meta = (mt_s0 + mt_s1 + mt_s2) / 3
    v2_meta = (v2[s0] + v2[s1] + v2[s2]) / 3
    ld_v2_s1 = np.log2(v2_s1 / v2_meta)
    ld_mt_s1 = np.log2(mt_s1 / mt_meta)
    cov_ld = float(np.mean(ld_v2_s1 * ld_mt_s1))
    print(f"  lam={lam}: Cov_raw={cov_raw:+.4e}  Cov_ld={cov_ld:+.4e}")

print("\ndone")
