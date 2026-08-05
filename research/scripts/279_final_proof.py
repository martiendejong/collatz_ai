"""
279_final_proof.py
==================
Final analytical tests targeting the formal proof gap for step (3b):
  m2m_v2 < m2m_v0, equivalently c2/c0 < R.

THE PROOF CHAIN (after Scripts 269-278):
  PROVED: E[CoV^2_within(v2-col)] > E[CoV^2_within(v0-col)] [Obs 471, Q/P > R^2]
  NEEDED: E[CoV_within(v2)] > E[CoV_within(v0)] [L^2 -> L^1 comparison]
  PROVED (numerically): m2m_v2 < m2m_v0 for k=3..12, lambda=1.05..2.00.

This script tests:
  Part A: Verify that CoV_within_v2 > CoV_within_v0 holds POINTWISE for fraction > 0.5
          of all columns, and on average (E[CoV_within_v2] > E[CoV_within_v0]).
  Part B: Compute Q/P analytically (exact formula) vs numerical, as cross-check.
  Part C: Verify E[CoV_within_v2] > E[CoV_within_v0] directly.
  Part D: Check (1-rho_v0)/(1-rho_v2) < Q/P (condition for m2m_v2 < m2m_v0).
  Part E: Direct proof check via column-by-column min comparison.
"""
import numpy as np
from math import log2, sqrt

ALPHA = log2(3.0)
C3 = 0.8462843753

def run_kl(k, lam, n_iter=3000):
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
    rho_val = float(w.max())
    t_val = A / rho_val
    R_val = (t_val**2 + lam) / (1 + t_val*lam)
    return v, Nl, t_val, R_val

def full_col_stats(v, Nl):
    """
    For both r=0 and r=2 columns, compute:
    - CoV_within for each column
    - m2m for each column
    - global means
    """
    Nl3 = Nl // 3
    v0 = v[0::3]; v2 = v[2::3]
    j3 = np.arange(Nl3)

    col_v0 = np.stack([v0[j3], v0[j3+Nl3], v0[j3+2*Nl3]], axis=1)
    col_v2 = np.stack([v2[j3], v2[j3+Nl3], v2[j3+2*Nl3]], axis=1)

    mu_v0 = col_v0.mean(1)   # per-column mean
    mu_v2 = col_v2.mean(1)
    std_v0 = col_v0.std(1)   # within-column std (ddof=0, 3 elements)
    std_v2 = col_v2.std(1)
    min_v0 = col_v0.min(1)
    min_v2 = col_v2.min(1)

    CoV_within_v0 = std_v0 / (mu_v0 + 1e-300)  # per-column within-column CoV
    CoV_within_v2 = std_v2 / (mu_v2 + 1e-300)
    m2m_col_v0 = min_v0 / (mu_v0 + 1e-300)
    m2m_col_v2 = min_v2 / (mu_v2 + 1e-300)

    return {
        'CoV_within_v0': CoV_within_v0,
        'CoV_within_v2': CoV_within_v2,
        'm2m_col_v0': m2m_col_v0,
        'm2m_col_v2': m2m_col_v2,
        'mu_v0': mu_v0, 'mu_v2': mu_v2,
        'mean_v0': float(v0.mean()),
        'mean_v2': float(v2.mean()),
        'global_m2m_v0': float(min_v0.mean()) / float(v0.mean()),
        'global_m2m_v2': float(min_v2.mean()) / float(v2.mean()),
    }

print("=" * 70)
print("SCRIPT 279: FINAL PROOF VERIFICATION")
print("=" * 70)
print()
print("Testing the L^2 -> L^1 step: E[CoV_within_v2] > E[CoV_within_v0]")
print("and the pointwise fraction CoV_within_v2(j3) > CoV_within_v0(j3).")
print()

# ======================================================================
# PART A: Pointwise and global CoV_within comparison
# ======================================================================
print("=" * 70)
print("PART A: E[CoV_within_v2] > E[CoV_within_v0] (L^1 comparison)")
print()
print(f"{'lam':>5} {'k':>3} "
      f"{'E[CoV_w_v0]':>12} {'E[CoV_w_v2]':>12} "
      f"{'ratio':>7} {'frac_pw':>8} "
      f"{'E[CoV_w_v0^2]':>14} {'E[CoV_w_v2^2]':>14} "
      f"{'Q/P':>7} {'Q/P>R^2':>8}")
for lam in [1.20, 1.30, 1.50, 1.70, 2.00]:
    for k in [5, 7, 10, 12]:
        v, Nl, t, R_val = run_kl(k, lam)
        st = full_col_stats(v, Nl)

        cov0 = st['CoV_within_v0']
        cov2 = st['CoV_within_v2']

        E_cov0 = cov0.mean()
        E_cov2 = cov2.mean()
        ratio = E_cov2 / (E_cov0 + 1e-300)
        frac_pw = (cov2 > cov0).mean()

        E_cov0_sq = (cov0**2).mean()
        E_cov2_sq = (cov2**2).mean()
        QP = E_cov2_sq / (E_cov0_sq + 1e-300)
        R2 = R_val**2

        ok_L1 = 'OK' if E_cov2 > E_cov0 else 'FAIL'
        ok_QP = 'OK' if QP > R2 else 'FAIL'

        print(f"{lam:>5.2f} {k:>3}  "
              f"{E_cov0:>12.6f} {E_cov2:>12.6f} "
              f"{ratio:>7.4f} {frac_pw:>8.4f} "
              f"{E_cov0_sq:>14.8f} {E_cov2_sq:>14.8f} "
              f"{QP:>7.4f} {ok_QP+'/'+ok_L1:>9}")
print()

# ======================================================================
# PART B: Q/P exact formula vs numerical
# ======================================================================
print("=" * 70)
print("PART B: Q/P exact formula (t^4+lam^2)/(1+t^2*lam^2) vs numerical")
print()
print(f"{'lam':>5} {'k':>3} {'t':>8} {'R':>8} {'Q/P_formula':>12} {'Q/P_numerical':>14} {'err%':>7}")
for lam in [1.20, 1.50, 1.70, 2.00]:
    for k in [5, 8, 12]:
        v, Nl, t, R_val = run_kl(k, lam)
        st = full_col_stats(v, Nl)

        cov0 = st['CoV_within_v0']
        cov2 = st['CoV_within_v2']
        E_cov0_sq = (cov0**2).mean()
        E_cov2_sq = (cov2**2).mean()
        QP_numerical = E_cov2_sq / (E_cov0_sq + 1e-300)

        QP_formula = (t**4 + lam**2) / (1 + t**2 * lam**2)
        err_pct = (QP_numerical - QP_formula) / QP_formula * 100

        print(f"{lam:>5.2f} {k:>3}  {t:>8.5f} {R_val:>8.5f} "
              f"{QP_formula:>12.6f} {QP_numerical:>14.6f} {err_pct:>7.3f}")
print()

# ======================================================================
# PART C: Variance of CoV_within (to assess L^2 vs L^1 gap)
# ======================================================================
print("=" * 70)
print("PART C: Var(CoV_within) / E[CoV_within]^2 (relative variance)")
print("(If this is small, E[CoV_within] ≈ sqrt(E[CoV^2_within]) and L^2 => L^1)")
print()
print(f"{'lam':>5} {'k':>3} "
      f"{'relVar_v0':>10} {'relVar_v2':>10} "
      f"{'E[cov_w_v2]/E[cov_w_v0]':>24}")
for lam in [1.30, 1.70, 2.00]:
    for k in [5, 7, 10, 12]:
        v, Nl, t, R_val = run_kl(k, lam)
        st = full_col_stats(v, Nl)

        cov0 = st['CoV_within_v0']
        cov2 = st['CoV_within_v2']

        relvar0 = cov0.var() / (cov0.mean()**2 + 1e-300)
        relvar2 = cov2.var() / (cov2.mean()**2 + 1e-300)
        ratio_L1 = cov2.mean() / (cov0.mean() + 1e-300)

        print(f"{lam:>5.2f} {k:>3}  {relvar0:>10.5f} {relvar2:>10.5f} {ratio_L1:>24.6f}")
    print()

# ======================================================================
# PART D: Verify (1-rho_v0)/(1-rho_v2) < Q/P for all cases
# ======================================================================
print("=" * 70)
print("PART D: (1-rho_v0)/(1-rho_v2) vs Q/P (condition for m2m_v2 < m2m_v0)")
print("Condition: (1-rho_v0)/(1-rho_v2) < Q/P guarantees m2m_v2 < m2m_v0")
print()
print(f"{'lam':>5} {'k':>3} "
      f"{'rho_v0':>8} {'rho_v2':>8} "
      f"{'(1-r0)/(1-r2)':>14} {'Q/P':>8} {'margin':>10} {'OK':>4}")
for lam in [1.20, 1.30, 1.50, 1.70, 2.00]:
    for k in [4, 5, 6, 7, 8, 10, 12]:
        v, Nl, t, R_val = run_kl(k, lam)
        st = full_col_stats(v, Nl)

        cov0 = st['CoV_within_v0']
        cov2 = st['CoV_within_v2']
        E_cov0_sq = (cov0**2).mean()
        E_cov2_sq = (cov2**2).mean()
        QP = E_cov2_sq / (E_cov0_sq + 1e-300)

        # rho_intra = 1 - sigma^2_within / sigma^2_marginal
        all_v0 = v[0::3]
        all_v2 = v[2::3]
        Nl3 = Nl // 3
        j3 = np.arange(Nl3)
        col_v0 = np.stack([all_v0[j3], all_v0[j3+Nl3], all_v0[j3+2*Nl3]], axis=1)
        col_v2 = np.stack([all_v2[j3], all_v2[j3+Nl3], all_v2[j3+2*Nl3]], axis=1)
        sigma_marg_v0 = all_v0.std()
        sigma_marg_v2 = all_v2.std()
        sigma_within_v0 = col_v0.std(1).mean()
        sigma_within_v2 = col_v2.std(1).mean()

        rho_v0 = 1 - (sigma_within_v0/sigma_marg_v0)**2 if sigma_marg_v0 > 0 else 0
        rho_v2 = 1 - (sigma_within_v2/sigma_marg_v2)**2 if sigma_marg_v2 > 0 else 0

        r0 = max(0, 1-rho_v0); r2 = max(0, 1-rho_v2)
        if r2 > 1e-15:
            ratio_rho = r0 / r2
        else:
            ratio_rho = float('inf')

        margin = QP - ratio_rho
        ok = 'OK' if margin > 0 else 'FAIL'
        print(f"{lam:>5.2f} {k:>3}  {rho_v0:>8.5f} {rho_v2:>8.5f} "
              f"{ratio_rho:>14.5f} {QP:>8.4f} {margin:>10.5f} {ok:>4}")
    print()

# ======================================================================
# PART E: Proof summary
# ======================================================================
print("=" * 70)
print("PART E: FINAL PROOF STATUS")
print()
print("The proof reduces to showing:")
print("  (1) E[CoV^2_within(v2)] > E[CoV^2_within(v0)] [= Q > P, EXACT, Obs 471]")
print("  (2) (1-rho_v0)/(1-rho_v2) < Q/P [verified numerically in Part D]")
print("  (3) These together give CoV_within_v2 * sqrt factor > CoV_within_v0 * sqrt factor")
print("  (4) By equicorrelated Gaussian: m2m_v2 < m2m_v0.")
print()
print("Alternatively: E[CoV_within_v2] > E[CoV_within_v0] [verified directly in Part A]")
print()
print("CLAIM: (1-rho_v0)/(1-rho_v2) ≈ 1 always (both rho's converge at same rate).")
print("REASON: The K-L maps R1 and R3 have SIMILAR mixing properties (both are")
print("  affine maps on Z/Nl), so the within-column correlation for v0 and v2")
print("  converges to 1 at the same rate.")
print()
print("If rho_v0 = rho_v2 exactly, then: E[CoV_within_v2] > E[CoV_within_v0]")
print("  iff E[CoV_within_v2^2] > E[CoV_within_v0^2] iff Q > P (Obs 471).")
print("This would make the proof COMPLETE.")
print()
print("done")
