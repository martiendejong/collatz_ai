"""
267_variance_proof.py
=====================
VARIANCE-BASED PROOF ATTEMPT for m2m_v2 <= m2m_v0.

INSIGHT: m2m = E[min(triplet)] / E[mean(triplet)].
For any positive random variable X, if Var(X) > 0 then E[min] < E[mean].
Specifically: E[min(X1,X2,X3)] / E[X] < 1 and decreases as Var(X)/E[X]^2 increases.

STRATEGY: Show within-column-triplet CoV (coefficient of variation) is larger for v2 than v0.
CoV_col_v2 > CoV_col_v0 => m2m_v2 < m2m_v0.

From K-L structure:
  v0[s] = (A/rho)*v2[sigma(s)] + (B1/rho)*cb[pi(s)]   [linear in v2, cb]
  v2[s] = (A^2/rho^2)*v0[phi(s)] + (B3/rho)*cb[R3(s)] [linear in v0, cb]

For v0 column triplet {s, s+Nl3, s+2Nl3} in sub-class r:
  - sigma permutes within r (maps r->r), so v2 inputs at sigma(sj) are in sub-class r
  - cb inputs at pi(sj) = (4sj)%Nl, class = sj%3 = r (pi preserves sub-class)
  => Both v2 and cb inputs vary across j within sub-class r

For v2 column triplet:
  - phi maps sc0->sc2, sc1->sc0, sc2->sc1 of v0 (phi is cross-class)
  - R3 maps sc0->sc1, sc1->sc0, sc2->sc2 of cb
  => Inputs come from DIFFERENT sub-classes across the triplet... NO:
  => Within a column triplet {s, s+Nl3, s+2Nl3} all have same sc (since Nl3/Nl = 1/3 and s%3=const)
  => phi(s)%3 = const for fixed s%3, so phi inputs are in ONE sub-class of v0

So BOTH v0 and v2 column-triplet inputs come from WITHIN their respective sub-classes.
The difference is T4 COEFFICIENT: A for v0, A^2/rho for v2 (smaller).

KEY VARIANCE DECOMPOSITION:
  Var_col(v0[s]) = (A/rho)^2 * Var_col(v2[sigma]) + (B1/rho)^2 * Var_col(cb[pi]) + cross
  Var_col(v2[s]) = (A^2/rho^2)^2 * Var_col(v0[phi]) + (B3/rho)^2 * Var_col(cb[R3]) + cross

CLAIM: Var_col(v2) > Var_col(v0) (adjusted for means).

This script verifies: CoV^2_col(v2) > CoV^2_col(v0) for all (k, lam).
"""
import numpy as np
from math import log2
import sys

ALPHA = log2(3.0)

def run_kl(k, lam, n_iter=None):
    if n_iter is None:
        n_iter = 700 + 100*max(0, k-8)
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

def analyze(k, lam):
    v, Nl, A, B1, B3 = run_kl(k, lam)
    Nl3 = Nl // 3
    j3 = np.arange(Nl3, dtype=np.int64)

    v0 = v[0::3]; v2 = v[2::3]

    # Column triplets for v0 and v2
    # triplet j: (v_type[j], v_type[j+Nl3], v_type[j+2*Nl3]) for j in [0, Nl3)
    col0 = np.stack([v0[j3], v0[j3+Nl3], v0[j3+2*Nl3]], axis=1)  # shape (Nl3, 3)
    col2 = np.stack([v2[j3], v2[j3+Nl3], v2[j3+2*Nl3]], axis=1)

    # Column means and mins
    mean_col0 = col0.mean(axis=1)
    mean_col2 = col2.mean(axis=1)
    min_col0  = col0.min(axis=1)
    min_col2  = col2.min(axis=1)

    # Column variance (within triplet)
    var_col0 = col0.var(axis=1)
    var_col2 = col2.var(axis=1)

    # CoV^2 = Var / Mean^2
    cov2_col0 = var_col0 / (mean_col0**2 + 1e-30)
    cov2_col2 = var_col2 / (mean_col2**2 + 1e-30)

    # Overall means
    mean_v0 = float(np.mean(v0))
    mean_v2 = float(np.mean(v2))
    c0 = float(np.mean(min_col0))
    c2 = float(np.mean(min_col2))
    m2m_v0 = c0 / mean_v0
    m2m_v2 = c2 / mean_v2

    # Mean CoV^2 within columns
    mean_cov2_v0 = float(np.mean(cov2_col0))
    mean_cov2_v2 = float(np.mean(cov2_col2))

    # DIRECT CHECK: min/mean ratio for each column
    mmr_col0 = min_col0 / (mean_col0 + 1e-30)
    mmr_col2 = min_col2 / (mean_col2 + 1e-30)
    mean_mmr_v0 = float(np.mean(mmr_col0))
    mean_mmr_v2 = float(np.mean(mmr_col2))

    # Variance of the within-column min/mean ratios
    # NOTE: m2m = E[min_col] / E[mean_col] (this is NOT E[min_col/mean_col] in general)
    # Confirm: E[min_col]/E[mean_col] = c0/mean_v0 (since E[mean_col] = mean_v0)
    # This equals m2m_v0 by definition.

    # Within-column CoV^2 > CoV^2 for v2?
    cov2_ok = mean_cov2_v2 > mean_cov2_v0

    # Per-sub-class breakdown
    sub = {}
    for sc in [0, 1, 2]:
        idx = j3[j3%3 == sc]
        if len(idx) == 0:
            continue
        c0_sc = np.stack([v0[idx], v0[idx+Nl3], v0[idx+2*Nl3]], axis=1)
        c2_sc = np.stack([v2[idx], v2[idx+Nl3], v2[idx+2*Nl3]], axis=1)
        cov2_v0_sc = float(np.mean(c0_sc.var(axis=1) / (c0_sc.mean(axis=1)**2 + 1e-30)))
        cov2_v2_sc = float(np.mean(c2_sc.var(axis=1) / (c2_sc.mean(axis=1)**2 + 1e-30)))
        sub[sc] = (cov2_v0_sc, cov2_v2_sc, cov2_v2_sc > cov2_v0_sc)

    return {
        'm2m_v0': m2m_v0, 'm2m_v2': m2m_v2, 'm2m_ok': m2m_v2 <= m2m_v0,
        'mean_cov2_v0': mean_cov2_v0, 'mean_cov2_v2': mean_cov2_v2,
        'cov2_ok': cov2_ok,
        'mean_mmr_v0': mean_mmr_v0, 'mean_mmr_v2': mean_mmr_v2,
        'sub': sub,
    }

print("267: Variance-based proof attempt for m2m_v2 <= m2m_v0")
print("="*70)

# Detailed at lam=1.70, k=10
lam, k = 1.70, 10
d = analyze(k, lam)
print(f"\nDetailed: lam={lam}, k={k}")
print(f"  m2m_v0 = {d['m2m_v0']:.6f}, m2m_v2 = {d['m2m_v2']:.6f}, OK: {d['m2m_ok']}")
print(f"  Mean within-col CoV^2 for v0: {d['mean_cov2_v0']:.6f}")
print(f"  Mean within-col CoV^2 for v2: {d['mean_cov2_v2']:.6f}")
print(f"  CoV^2(v2) > CoV^2(v0): {d['cov2_ok']}")
print(f"  Mean within-col min/mean for v0: {d['mean_mmr_v0']:.6f}")
print(f"  Mean within-col min/mean for v2: {d['mean_mmr_v2']:.6f}")
print()
for sc, (cv0, cv2, ok) in d['sub'].items():
    print(f"  sc={sc}: CoV^2_v0={cv0:.5f}, CoV^2_v2={cv2:.5f}, v2>v0: {ok}")

print(f"\nScan (k=8,10): CoV^2(v2) > CoV^2(v0) for all lambda?")
print(f"{'lam':>6} {'k':>3} {'CoV2_v0':>10} {'CoV2_v2':>10} {'v2>v0':>8} {'m2m_ok':>8}")
for lam in [1.30, 1.50, 1.70, 1.90, 2.00]:
    for k in [8, 10]:
        d = analyze(k, lam)
        print(f"lam={lam:.2f} k={k:>2} {d['mean_cov2_v0']:>10.6f} {d['mean_cov2_v2']:>10.6f} {str(d['cov2_ok']):>8} {str(d['m2m_ok']):>8}")
    sys.stdout.flush()

print("""
ANALYTICAL ARGUMENT (if CoV^2 ordering is verified):

For any positive column-triplet (X1, X2, X3) with mean mu and CoV^2 = sigma^2/mu^2:
  E[min] / E[mean] = m2m = 1 - f(CoV^2) approximately, where f is increasing.

Exact relation (for ordered stats of 3 values):
  E[min(U,V,W)] / E[mean(U,V,W)] = 1 - Cov(indicator, value) / E[mean]

If CoV^2(v2 col) > CoV^2(v0 col) for ALL columns simultaneously:
  => min/mean is lower for v2 columns than v0 columns (more spread => lower min)
  => E[min_v2_col] / E[mean_v2_col] < E[min_v0_col] / E[mean_v0_col]
  => m2m_v2 < m2m_v0. QED.

Note: m2m = E[min_col] / E[mean_col] = E[min_col] / mean_v
And E[min_col] = E[min(X1,X2,X3)] where (X1,X2,X3) is a column-triplet.

For an EXACT bound: for any positive (X1,X2,X3):
  min(X1,X2,X3) / mean(X1,X2,X3) <= 1 with equality iff X1=X2=X3.
  Moreover, min/mean >= 1 - sqrt(2/3) * sqrt(Var/Mean^2) approximately
  (since the min of 3 equidistributed normals is ~1 - 0.85*sigma/mu).

If we can prove CoV^2(v2) > CoV^2(v0) for each column, the ordering follows.
""")
print("done")
