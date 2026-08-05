"""
265_pointwise_order.py
======================
ATTEMPT: Prove m2m_v2 <= m2m_v0 analytically.

APPROACH: Study whether there's a POINTWISE ratio bound v2[s]/v0[sigma(s)]
that creates the ordering.

Key: c2/c0 = (rho*a2_v0 - A*a2_v2)/(rho*a0_v0 - A*a0_v2) = B1*c2/(B1*c0) -- circular!
     (from y_r = (A/rho)*x_r + (B1/rho)*c_r => rho*y_r - A*x_r = B1*c_r)

So the gap c2/c0 <= R cannot be derived from the 6-variable MEAN equations alone --
it requires distributional information (variance, not just means).

ALTERNATIVE: Does the pointwise bound v2[s] <= K*v0[pi_inv(s)] hold for some K < R?
If v2[s]/v0[s'] <= R for all s, s' on same T4-chain, then column-mins also bound by R.

Check 1: Ratio v2[s]/v0[s'] for pairs (s,s') linked via T4 chain.
Check 2: Distribution of v2[s]/v0[s] for all s.
Check 3: Does max(v2)/max(v0) == 1 (by normalization) while min(v2)/min(v0) < R?
Check 4: Within each sub-class, is E[v2]/E[v0] == R? (Should be by definition.)
Check 5: Is the WITHIN-SUBCLASS ratio c2_subclass/a2_v2 <= c0_subclass/a0_v0?
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
    s = np.arange(Nl, dtype=np.int64)

    v0 = v[0::3]  # r=0 nodes, indexed by s in [0,Nl)
    v2 = v[2::3]  # r=2 nodes, indexed by s in [0,Nl)

    mean_v0 = float(np.mean(v0))
    mean_v2 = float(np.mean(v2))
    R = mean_v2/mean_v0

    # rho
    v1 = v[1::3]
    sigma1 = (4*s+2) % Nl
    rho = A / float(np.mean(v1 / v0[sigma1]))

    # Per-sub-class analysis
    Nl3 = Nl//3
    j3 = np.arange(Nl3, dtype=np.int64)

    results = {}
    for sc in [0, 1, 2]:
        mask = (s % 3 == sc)
        v0_sc = v0[mask]
        v2_sc = v2[mask]
        mean_v0_sc = float(np.mean(v0_sc))
        mean_v2_sc = float(np.mean(v2_sc))
        R_sc = mean_v2_sc/mean_v0_sc if mean_v0_sc > 0 else float('inf')

        # Column-min for this sub-class
        # j values for this sub-class in column-min: j in [0,Nl3) with j%3==sc
        j_sc = j3[j3%3==sc]  # sub-class sc indices within first Nl3
        cb_v0_col_sc = np.minimum(np.minimum(v0[j_sc], v0[j_sc+Nl3]), v0[j_sc+2*Nl3])
        cb_v2_col_sc = np.minimum(np.minimum(v2[j_sc], v2[j_sc+Nl3]), v2[j_sc+2*Nl3])
        c0_sc = float(np.mean(cb_v0_col_sc))
        c2_sc = float(np.mean(cb_v2_col_sc))
        m2m_v0_sc = c0_sc/mean_v0_sc if mean_v0_sc > 0 else float('nan')
        m2m_v2_sc = c2_sc/mean_v2_sc if mean_v2_sc > 0 else float('nan')
        R_colmin = c2_sc/c0_sc if c0_sc > 0 else float('inf')

        results[sc] = {
            'mean_v0': mean_v0_sc, 'mean_v2': mean_v2_sc, 'R_mean': R_sc,
            'c0': c0_sc, 'c2': c2_sc, 'R_colmin': R_colmin,
            'm2m_v0': m2m_v0_sc, 'm2m_v2': m2m_v2_sc,
            'm2m_v2<=m2m_v0': m2m_v2_sc <= m2m_v0_sc,
            'R_colmin<=R_mean': R_colmin <= R_sc,
        }

    # Overall
    Nl3 = Nl//3
    cb_v0_col = np.minimum(np.minimum(v0[j3], v0[j3+Nl3]), v0[j3+2*Nl3])
    cb_v2_col = np.minimum(np.minimum(v2[j3], v2[j3+Nl3]), v2[j3+2*Nl3])
    c0 = float(np.mean(cb_v0_col))
    c2 = float(np.mean(cb_v2_col))
    m2m_v0 = c0/mean_v0
    m2m_v2 = c2/mean_v2

    return {
        'R': R, 'c0': c0, 'c2': c2, 'c2_c0': c2/c0,
        'm2m_v0': m2m_v0, 'm2m_v2': m2m_v2,
        'm2m_ok': m2m_v2 <= m2m_v0,
        'subclass': results,
    }

print("265: Per-sub-class m2m analysis")
print("="*75)

lam, k = 1.70, 10
d = analyze(k, lam)
print(f"\nlam={lam}, k={k}")
print(f"  Overall: c2/c0={d['c2_c0']:.5f}, R={d['R']:.5f}, m2m_v2={d['m2m_v2']:.5f}, m2m_v0={d['m2m_v0']:.5f}")
print(f"  m2m_v2 <= m2m_v0: {d['m2m_ok']}")
print()
print(f"  {'sc':>4}  {'mean_v0':>10}  {'mean_v2':>10}  {'R_mean':>8}  {'c0_sc':>10}  {'c2_sc':>10}  {'R_col':>8}  {'m2m_v0':>8}  {'m2m_v2':>8}  {'m2m_ok':>8}  {'R_col<=R':>10}")
for sc in [0, 1, 2]:
    r = d['subclass'][sc]
    print(f"  sc={sc}  {r['mean_v0']:>10.5f}  {r['mean_v2']:>10.5f}  {r['R_mean']:>8.5f}  {r['c0']:>10.5f}  {r['c2']:>10.5f}  {r['R_colmin']:>8.5f}  {r['m2m_v0']:>8.5f}  {r['m2m_v2']:>8.5f}  {str(r['m2m_v2<=m2m_v0']):>8}  {str(r['R_colmin<=R_mean']):>10}")

print()

# Lambda scan k=8
print(f"\nLambda scan k=8: per-sub-class R_colmin vs R_mean")
print(f"{'lam':>6}  {'sc0: R_col<=R':>14}  {'sc1: R_col<=R':>14}  {'sc2: R_col<=R':>14}  {'overall m2m_ok':>16}")
k = 8
for lam in [1.30, 1.50, 1.70, 1.90, 2.00]:
    d = analyze(k, lam)
    sc0_ok = d['subclass'][0]['R_colmin<=R_mean']
    sc1_ok = d['subclass'][1]['R_colmin<=R_mean']
    sc2_ok = d['subclass'][2]['R_colmin<=R_mean']
    print(f"lam={lam:.2f}  {str(sc0_ok):>14}  {str(sc1_ok):>14}  {str(sc2_ok):>14}  {str(d['m2m_ok']):>16}")
    sys.stdout.flush()

print(f"""
KEY OBSERVATION:
  m2m inequality m2m_v2 <= m2m_v0 holds OVERALL but may fail PER-SUB-CLASS.
  If it holds per-sub-class AND sub-class means of v2 >= corresponding v0 means,
  then the overall inequality follows by weighted averaging.

  If NOT per-sub-class: the inequality is a MIXING EFFECT across sub-classes.
  In that case: need a covariance/rearrangement argument.
""")

print("done")
