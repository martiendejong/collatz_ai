"""
277_deeper_analysis.py
======================
Deeper analytical investigation of step (3b): c2/c0 < R.

MAIN INSIGHT UNDER INVESTIGATION:
  The sub-additive inequality gives c0 >= (A/rho)*c2 + (B1/rho)*cc0'.
  This implies c2/c0 <= 1/t - lam^alpha * cc0'/c0.
  For this to give c2/c0 <= R, need cc0'/c0 >= mean_cb/mean_v0.

  However: cc0' -> 0 as k -> inf (global minimum of 9-element blocks -> 0).
  So the sub-additive bound WEAKENS to c2/c0 <= 1/t (not useful for proving <= R).

  THEREFORE: the sub-additive bound alone cannot prove c2/c0 < R.

ALTERNATIVE APPROACH:
  The correct proof must use the SPECIFIC STRUCTURE of the K-L eigenvector,
  not just general sub-additivity.

  This script investigates:
  (A) The EXACT column-by-column contribution to c2/c0 vs R.
  (B) A "BALANCE EQUATION" relating c2/c0 across depths.
  (C) The SECOND-MOMENT bound for triplet minima.
  (D) Exact computation of c2/c0 for k=3 (N=9) symbolically.
"""
import numpy as np
from math import log2, sqrt
import sympy as sp
from fractions import Fraction

ALPHA = log2(3.0)

def run_kl(k, lam, n_iter=6000):
    A  = lam**-2.0; B1 = lam**(ALPHA-2.0); B3 = lam**(ALPHA-1.0)
    N  = 3**(k-1); Nl = N//3
    i  = np.arange(N, dtype=np.int64)
    T4 = (4*i+2) % N
    s_arr, r_arr = np.divmod(i, 3)
    m0_mask, m2_mask = (r_arr==0), (r_arr==2)
    R1 = (4*s_arr) % Nl; R3 = (2*s_arr+1) % Nl
    v = np.ones(N, dtype=np.float64)
    for it in range(n_iter):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        w  = A*v[T4]
        w[m2_mask] += B3*cb[R3[m2_mask]]
        w[m0_mask] += B1*cb[R1[m0_mask]]
        vmax = w.max()
        w /= vmax
        if it > 200 and it % 1000 == 0:
            if np.abs(w-v).max() < 1e-14:
                break
        v = w
    return v, Nl, A, B1, B3, T4, s_arr, r_arr, R1, R3, m0_mask, m2_mask

# ===========================================================================
# PART A: PER-COLUMN ANALYSIS
# Compare at column-level: (min(v2-col j3)/mean(v2-col j3)) vs
#                          (min(v0-col j3)/mean(v0-col j3))
# and how these WEIGHT the global c2/c0 vs R.
# ===========================================================================
print("="*72)
print("PART A: Per-column analysis (column-level contributions)")
print("="*72)
print()

for lam in [1.50, 1.70, 2.00]:
    for k in [6, 10]:
        v, Nl, A, B1, B3, T4, s_arr, r_arr, R1, R3, m0_mask, m2_mask = run_kl(k, lam)
        N = 3**(k-1)
        Nl3 = Nl // 3
        v0 = v[0::3]; v2 = v[2::3]
        j3 = np.arange(Nl3)

        col_v0 = np.stack([v0[j3], v0[j3+Nl3], v0[j3+2*Nl3]], axis=1)
        col_v2 = np.stack([v2[j3], v2[j3+Nl3], v2[j3+2*Nl3]], axis=1)

        min_v0 = col_v0.min(1); mean_v0_col = col_v0.mean(1)
        min_v2 = col_v2.min(1); mean_v2_col = col_v2.mean(1)

        c0 = min_v0.mean(); c2 = min_v2.mean()
        mv0 = v0.mean(); mv2 = v2.mean()
        R_val = mv2 / mv0

        # Column-level min/mean ratios
        m2m_v0_col = min_v0 / mean_v0_col  # per-column m2m for v0
        m2m_v2_col = min_v2 / mean_v2_col  # per-column m2m for v2

        # Global m2m ratios
        m2m_v0_global = c0 / mv0
        m2m_v2_global = c2 / mv2

        # Covariance decomposition:
        # E[min(v0-col)] = E[m2m_v0_col * mean_v0_col]
        # = E[m2m_v0_col] * E[mean_v0_col] + cov(m2m_v0_col, mean_v0_col)
        # E[mean_v0_col] = mv0
        cov_m2m_mean_v0 = np.cov(m2m_v0_col, mean_v0_col)[0,1]
        cov_m2m_mean_v2 = np.cov(m2m_v2_col, mean_v2_col)[0,1]

        # c2/c0 decomposition:
        # c2/c0 = E[m2m_v2_col * mean_v2_col] / E[m2m_v0_col * mean_v0_col]
        # R = mv2/mv0 = E[mean_v2_col] / E[mean_v0_col]

        # Weighted average: c2/c0 = sum(min_v2) / sum(min_v0)
        # = [sum(m2m_v2_col * mean_v2_col)] / [sum(m2m_v0_col * mean_v0_col)]

        # If m2m_v2_col/m2m_v0_col = const = q for all j3:
        # then c2/c0 = q * sum(mean_v2_col) / sum(mean_v0_col) = q * R.
        # So c2/c0 < R iff q < 1 iff m2m_v2_col < m2m_v0_col (on average, weighted).

        # Weighted mean of m2m_v2_col with weights mean_v2_col:
        w2 = mean_v2_col / mean_v2_col.sum()
        w0 = mean_v0_col / mean_v0_col.sum()
        wm_v2 = (m2m_v2_col * w2).sum()  # weighted m2m for v2
        wm_v0 = (m2m_v0_col * w0).sum()  # weighted m2m for v0

        print(f"lam={lam:.2f}, k={k}:")
        print(f"  c2/c0={c2/c0:.8f}, R={R_val:.8f}, diff={R_val-c2/c0:.2e}")
        print(f"  m2m_v0(global)={m2m_v0_global:.8f}, m2m_v2(global)={m2m_v2_global:.8f}")
        print(f"  weighted_m2m_v0={wm_v0:.8f}, weighted_m2m_v2={wm_v2:.8f}")
        print(f"  c2/c0 = wm_v2 * R?  wm_v2*R={wm_v2*R_val:.8f} vs c2/c0={c2/c0:.8f}")
        print(f"  c2/c0 = wm_v2/wm_v0 * R?  {(wm_v2/wm_v0)*R_val:.8f}")
        print(f"  => c2/c0 < R iff wm_v2 < wm_v0 (weighted m2m comparison)")
        print(f"  wm_v2 < wm_v0: {wm_v2 < wm_v0}")
        print(f"  cov(m2m_v0_col, mean_v0_col)={cov_m2m_mean_v0:.6e}")
        print(f"  cov(m2m_v2_col, mean_v2_col)={cov_m2m_mean_v2:.6e}")
        print()

# ===========================================================================
# PART B: The "balance equation" approach
# From the K-L equation: c2/c0 = R * (something).
# Derive what that "something" is analytically.
# ===========================================================================
print("="*72)
print("PART B: Analytical expression for c2/c0 via K-L balance equation")
print("="*72)
print("""
From K-L: v0[i] = t * (v2[T4(i)] + lam^alpha * cb[R1(i)])
         v2[i] = t * (t*v0[T'(i)] + lam^{alpha+1} * cb[R3(i)])

Taking column minima and means, using sub-additivity:
  c0 >= t * (c2 + lam^alpha * cc0')
  c2 >= t * (t*c0 + lam^{alpha+1} * cc2')

where cc0' = E[min(cb-col_for_v0)] and cc2' = E[min(cb-col_for_v2)].

The SUB-ADDITIVE BOUNDS give c2/c0 <= 1/t (too weak).

THE CORRECT IDENTITY (not inequality) would be:
  c0 = t * (c2 + lam^alpha * cc0') + SLACK0
  c2 = t * (t*c0 + lam^{alpha+1} * cc2') + SLACK2

where SLACK0, SLACK2 >= 0 measure the non-alignment of minimizers.

FROM THE BALANCE: c2/c0 = R - (correction related to SLACK terms).
Since SLACK >= 0, this gives c2/c0 <= R when the correction is non-negative.

QUESTION: Is the correction always non-negative?

The correction involves: how much the minimizer of v2[T4-col]
DIFFERS from the minimizer of cb[R1-col] within the same v0-column.

When they COINCIDE: SLACK0 = max possible (equality holds in sub-add).
When they DIFFER: SLACK0 = 0 (no slack -- but this can't happen at finite k
since the T4 permutation mixes v2 indices from different scales).

NUMERICAL VERIFICATION of SLACK positivity:
""")

for lam in [1.50, 1.70, 2.00]:
    for k in [5, 8, 10]:
        v, Nl, A, B1, B3, T4, s_arr, r_arr, R1, R3, m0_mask, m2_mask = run_kl(k, lam)
        N = 3**(k-1); Nl3 = Nl//3
        t = A / (v.max())  # approximate rho = v.max()=1, so t = A
        # Actually rho*v = F(v), and v is normalized to max=1.
        # Let's compute rho properly.
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        w = A*v[T4]
        w[m2_mask] += B3*cb[R3[m2_mask]]
        w[m0_mask] += B1*cb[R1[m0_mask]]
        rho = w.max()  # since v normalized to max=1
        t_val = A/rho

        v0 = v[0::3]; v2 = v[2::3]
        j3 = np.arange(Nl3)
        col_v0 = np.stack([v0[j3], v0[j3+Nl3], v0[j3+2*Nl3]], axis=1)
        col_v2 = np.stack([v2[j3], v2[j3+Nl3], v2[j3+2*Nl3]], axis=1)
        c0 = float(col_v0.min(1).mean()); c2 = float(col_v2.min(1).mean())
        R_val = v2.mean()/v0.mean()

        # Compute the "slack" for v0 columns:
        # slack[j3] = min_s(v0[j3+s*Nl3]) - (A/rho)*min_s(v2[T4(j3+s*Nl3)])
        #           - (B1/rho)*min_s(cb[R1(j3+s*Nl3)])
        lam_alpha = lam**ALPHA

        # v2[T4-column at j3]: T4 maps v0-positions to v2-positions
        # T4(j3+s*Nl3) for s=0,1,2 gives three v2 positions
        T4_col = np.stack([T4[3*(j3+s*Nl3)] for s in range(3)], axis=1)  # shape (Nl3, 3)
        # Wait -- T4 is defined on ALL indices 0..N-1.
        # For r=0 elements: i = 3*(j3+s*Nl3), so T4[i] = T4[3*(j3+s*Nl3)].
        # v2[T4[i]] where i is an r=0 element.
        # v2 = v[2::3], so v2[j] = v[3j+2].
        # v2[T4[i]] = v[T4[i]] (since T4[i] has T4[i]%3=2 for r=0 inputs).

        # Let me compute directly:
        v0_positions = 3*(j3[:,None] + np.array([0,1,2])*Nl3)  # shape (Nl3, 3)
        v2_T4_vals = v[T4[v0_positions]]  # v2 at T4 of v0-positions
        cb_R1_vals = cb[R1[j3[:,None] + np.array([0,1,2])*Nl3]]  # cb at R1 of s-indices

        min_v0_col = col_v0.min(1)
        min_v2T4_col = v2_T4_vals.min(1)
        min_cbR1_col = cb_R1_vals.min(1)

        slack0 = min_v0_col - (A/rho)*min_v2T4_col - (B1/rho)*min_cbR1_col
        slack0_mean = slack0.mean()
        slack0_min = slack0.min()
        slack0_frac_positive = (slack0 > 0).mean()

        print(f"lam={lam:.2f}, k={k}: SLACK0 mean={slack0_mean:.6e}, "
              f"min={slack0_min:.2e}, frac>0={slack0_frac_positive:.3f}")

print()

# ===========================================================================
# PART C: Second-moment bound for triplet minimum
# For independent X1,X2,X3 with mean mu and variance sigma^2:
# E[min(X1,X2,X3)] >= mu - sqrt(2/3) * sigma (classical bound).
# Higher variance => lower minimum.
# Since Var(v2-col) > Var(v0-col) (Obs 471/476):
# E[min(v2-col)] / E[v2-col-mean] <= E[min(v0-col)] / E[v0-col-mean]
# (more variance => lower min/mean ratio).
# This is EXACTLY the direction we need!
# BUT: the bound requires INDEPENDENCE, not satisfied here.
# ===========================================================================
print("="*72)
print("PART C: Second-moment bound and variance-minimum relationship")
print("="*72)
print("""
For a triplet (X1,X2,X3) with mean mu and variance sigma^2 (per element):
  E[min] = mu - sqrt(2/3) * sigma + O(skewness terms)   [for symmetric dist]
  m2m = E[min]/mu = 1 - (1/mu) * sqrt(2/3) * sigma + ...
       = 1 - sqrt(2/3) * CoV + ...   [CoV = sigma/mu]

Since CoV(v2-col) > CoV(v0-col) (from Obs 471):
  m2m(v2) approx 1 - sqrt(2/3)*CoV(v2) < 1 - sqrt(2/3)*CoV(v0) approx m2m(v0).

This argument is APPROXIMATE (holds for symmetric distributions, not general).
Let's verify how well this approximation works:
""")

for lam in [1.50, 1.70, 2.00]:
    k = 10
    v, Nl, A, B1, B3, T4, s_arr, r_arr, R1, R3, m0_mask, m2_mask = run_kl(k, lam)
    N = 3**(k-1); Nl3 = Nl//3
    v0 = v[0::3]; v2 = v[2::3]
    j3 = np.arange(Nl3)
    col_v0 = np.stack([v0[j3], v0[j3+Nl3], v0[j3+2*Nl3]], axis=1)
    col_v2 = np.stack([v2[j3], v2[j3+Nl3], v2[j3+2*Nl3]], axis=1)

    # Per-column CoV for v0 and v2
    cov_v0_col = col_v0.std(1) / col_v0.mean(1)  # per-column CoV
    cov_v2_col = col_v2.std(1) / col_v2.mean(1)

    # Predicted m2m from second-moment bound:
    pred_m2m_v0 = 1 - np.sqrt(2/3) * cov_v0_col.mean()
    pred_m2m_v2 = 1 - np.sqrt(2/3) * cov_v2_col.mean()

    # Actual m2m:
    act_m2m_v0 = col_v0.min(1).mean() / col_v0.mean(1).mean()
    act_m2m_v2 = col_v2.min(1).mean() / col_v2.mean(1).mean()

    R_val = v2.mean()/v0.mean()
    c2_c0 = col_v2.min(1).mean() / col_v0.min(1).mean()

    print(f"lam={lam:.2f}, k={k}:")
    print(f"  Predicted m2m_v0={pred_m2m_v0:.6f}, actual={act_m2m_v0:.6f}")
    print(f"  Predicted m2m_v2={pred_m2m_v2:.6f}, actual={act_m2m_v2:.6f}")
    print(f"  Predicted: m2m_v2 < m2m_v0? {pred_m2m_v2 < pred_m2m_v0}")
    print(f"  Actual:    m2m_v2 < m2m_v0? {act_m2m_v2 < act_m2m_v0}")
    print(f"  c2/c0={c2_c0:.8f}, R={R_val:.8f}")
    print()

# ===========================================================================
# PART D: Direct symbolic computation for k=3 (N=9)
# Prove c2/c0 < R analytically for k=3 by computing eigenvector symbolically.
# ===========================================================================
print("="*72)
print("PART D: Exact computation k=3 (N=9) for lambda=2")
print("="*72)
print()

# For lambda=2, k=3:
# N=9, Nl=3, T4=(4i+2)%9, m0=(i%3==0), m2=(i%3==2)
# s_arr = i//3, R1=(4*s)%3, R3=(2*s+1)%3

lam = 2.0
A = lam**-2.0; B1 = lam**(ALPHA-2.0); B3 = lam**(ALPHA-1.0)
print(f"A={A:.6f}, B1={B1:.6f}, B3={B3:.6f}")
print(f"B1 = 3/4 = {3/4:.6f}, B3 = 3/2 = {3/2:.6f}")

# Compute numerically at high precision for k=3
v, Nl, A, B1, B3, T4, s_arr, r_arr, R1, R3, m0_mask, m2_mask = run_kl(3, lam, n_iter=50000)
N = 9; Nl = 3; Nl3 = 1

print(f"\nNumerical Perron eigenvector (k=3, lam=2, N=9, normalized to max=1):")
for i in range(N):
    cb_i = np.minimum(np.minimum(v[0], v[3]), v[6])  # Nl=3, so cb[j]=min(v[j],v[j+3],v[j+6])
    r_i = i % 3; s_i = i // 3
    print(f"  v[{i}] = {v[i]:.12f}  (r={r_i}, s={s_i}, T4={T4[i]})")

# Compute cb
cb3 = np.minimum(np.minimum(v[:3], v[3:6]), v[6:])
print(f"\ncb = {cb3}")

# Compute c2/c0 and R
v0_k3 = v[0::3]; v2_k3 = v[2::3]
# With Nl3=1, there's only 1 column triplet:
col_v0_k3 = np.array([[v0_k3[0], v0_k3[0], v0_k3[0]]])  # trivially 1x3
# Wait: v0 has Nl=3 elements, Nl3=Nl//3=1. So j3=[0].
# col_v0[j3=0] = {v0[0], v0[1], v0[2]} = {v[0], v[3], v[6]}
col_v0_k3 = np.array([[v0_k3[0], v0_k3[1], v0_k3[2]]])  # 1 column triplet
col_v2_k3 = np.array([[v2_k3[0], v2_k3[1], v2_k3[2]]])

c0_k3 = float(col_v0_k3.min(1).mean())
c2_k3 = float(col_v2_k3.min(1).mean())
R_k3 = float(v2_k3.mean() / v0_k3.mean())
print(f"\nc0 = min(v[0],v[3],v[6]) = {c0_k3:.12f}")
print(f"c2 = min(v[2],v[5],v[8]) = {c2_k3:.12f}")
print(f"R = mean(v2)/mean(v0) = {R_k3:.12f}")
print(f"c2/c0 = {c2_k3/c0_k3:.12f}")
print(f"c2/c0 < R: {c2_k3/c0_k3 < R_k3}")
print(f"Margin R - c2/c0 = {R_k3 - c2_k3/c0_k3:.2e}")

# The exact eigenvector for k=3 satisfies a polynomial of degree 9.
# Let's verify the K-L fixed point equation at these values.
rho_est = (A*v[T4] + B3*cb3[R3]*m2_mask.astype(float) + B1*cb3[R1]*m0_mask.astype(float)).max()
print(f"\nEstimated rho = {rho_est:.12f}")
print(f"t = A/rho = {A/rho_est:.12f}")
print(f"R formula (t^2+lam)/(1+t*lam) = {(A/rho_est)**2+lam:.6f}/{1+(A/rho_est)*lam:.6f}")
R_formula = ((A/rho_est)**2 + lam) / (1 + (A/rho_est)*lam)
print(f"  = {R_formula:.12f}")
print(f"R from eigenvector = {R_k3:.12f}")
print(f"Difference = {abs(R_formula-R_k3):.2e}")

# ===========================================================================
# PART E: THE KEY OBSERVATION
# ===========================================================================
print()
print("="*72)
print("PART E: Key theoretical insight (Obs 479)")
print("="*72)
print("""
FUNDAMENTAL STRUCTURE OF c2/c0 < R:

The K-L eigenvector at depth k satisfies:
  v0-col triplet: {(A/rho)*v2[T4_0] + (B1/rho)*cb[R1_0],
                   (A/rho)*v2[T4_1] + (B1/rho)*cb[R1_1],
                   (A/rho)*v2[T4_2] + (B1/rho)*cb[R1_2]}

  v2-col triplet: {(A^2/rho^2)*v0[T'_0] + (B3/rho)*cb[R3_0],
                   (A^2/rho^2)*v0[T'_1] + (B3/rho)*cb[R3_1],
                   (A^2/rho^2)*v0[T'_2] + (B3/rho)*cb[R3_2]}

The CRUCIAL ASYMMETRY:
  * v0-col has coefficient ratio (B1/rho) / (A/rho) = B1/A = lam^alpha.
  * v2-col has coefficient ratio (B3/rho) / (A^2/rho^2) = B3*rho/A^2 = lam^(alpha+1)*rho/A
    = lam^(alpha+1)/t.

So v2-col has a PROPORTIONALLY LARGER cb contribution (by factor lam/t > lam > 1).

This means: the cb-induced variation is lam/t times larger in v2-col than in v0-col.
The cb-variance "overwhelms" the backbone variance for v2.

From Obs 471: Q/P = (t^4+lam^2)/(1+t^2*lam^2) > R^2 = ((t^2+lam)/(1+t*lam))^2.
This is the ANALYTICAL PROOF that V(v2-col) > V(v0-col) when expressed in terms
of the Perron eigenvector coefficients.

THE VARIANCE DIFFERENTIAL DIRECTLY PREDICTS m2m_v2 < m2m_v0:
  Higher within-column variance => lower min/mean ratio.
  Since Var(v2-col) > Var(v0-col) GLOBALLY (Obs 471), we expect m2m_v2 < m2m_v0.

THE GAP IN THE PROOF:
  The "higher variance => lower m2m" connection is NOT universally rigorous.
  It holds for:
  (a) log-normal distributions: FAILS (error too large, Script 271).
  (b) Second-order stochastic dominance: REQUIRES more than just variance comparison.
  (c) The specific K-L joint distribution: needs detailed analysis.

WHAT WOULD CLOSE THE PROOF:
  A function f: R -> R such that m2m(col) = f(Var(col)) + (lower-order terms),
  with f strictly decreasing. This would convert the analytical variance bound
  (Obs 471) into an analytical m2m bound.

  For GAUSSIAN distributions: E[min(X1,X2,X3)] / E[Xi] = mu - sqrt(2/3)*sigma,
  so m2m = 1 - sqrt(2/3)*CoV. The functional IS f(Var) = 1 - sqrt(2/3)*sqrt(Var/mu^2).
  This gives m2m_v2 < m2m_v0 DIRECTLY from Var(v2) > Var(v0).

  QUESTION: Are the K-L column distributions close enough to Gaussian for this to hold?
""")

# Check: how Gaussian are the column distributions?
print("Checking Gaussianity of K-L column distributions:")
from scipy import stats as scipy_stats

for lam in [1.50, 2.00]:
    k = 10
    v, Nl, A, B1, B3, T4, s_arr, r_arr, R1, R3, m0_mask, m2_mask = run_kl(k, lam)
    N = 3**(k-1); Nl3 = Nl//3
    v0 = v[0::3]; v2 = v[2::3]
    j3 = np.arange(Nl3)
    col_v0 = np.stack([v0[j3], v0[j3+Nl3], v0[j3+2*Nl3]], axis=1)
    col_v2 = np.stack([v2[j3], v2[j3+Nl3], v2[j3+2*Nl3]], axis=1)

    # For each column, compute the actual min and the predicted min (Gaussian approximation)
    for name, col in [('v0', col_v0), ('v2', col_v2)]:
        mu = col.mean(1); sigma = col.std(1)
        # Predicted min for Gaussian: mu - sqrt(2/3)*sigma (but this is for 3 iid)
        # For order statistics: E[min of 3 iid N(0,1)] = -sqrt(3/2) = -sqrt(6)/2
        # So E[min] = mu - sqrt(6)/2 * sigma/sqrt(1) ... check
        # For 3 iid N(mu, sigma^2): E[min] = mu - sqrt(6)/(2) * sigma? Actually:
        # E[X(1:3)] where X(1:3) is the minimum of 3 iid N(0,1) = -3*phi(0)/something
        # Standard result: E[min(X1,X2,X3)] = mu - sigma * E[min(Z1,Z2,Z3)] where Zi~N(0,1)
        # E[min(Z1,Z2,Z3)] = -3*phi(0)*sqrt(2) ... actually:
        # For n=3 iid N(0,1): E[Z(1:3)] = -3*sqrt(3)/(2*sqrt(pi)) = -3*sqrt(3)/(2*sqrt(pi))
        E_min_standard = -3*np.sqrt(3) / (2*np.sqrt(np.pi))  # = -1.4142... /sqrt(pi) * ...
        # Actually the correct value is:
        # E[min of 3 iid N(0,1)] = integral from -inf to inf of x * 3*(1-Phi(x))^2 * phi(x) dx
        # = -sqrt(3)/sqrt(pi) * (some constant)
        # Let me just use the numerical value:
        from scipy.stats import norm
        def E_min_n(n):
            # E[min of n iid N(0,1)] via numerical integration
            from scipy.integrate import quad
            def integrand(x):
                return x * n * (1 - norm.cdf(x))**(n-1) * norm.pdf(x)
            result, _ = quad(integrand, -10, 10)
            return result
        c_n3 = E_min_n(3)  # should be about -0.8463

        pred_min = mu + c_n3 * sigma  # E[min] predicted by Gaussian
        act_min = col.min(1)
        resid = act_min - pred_min

        corr = np.corrcoef(pred_min, act_min)[0,1]
        rmse = np.sqrt((resid**2).mean())
        rel_rmse = rmse / mu.mean()

        print(f"  lam={lam:.2f}, k={k}, {name}: corr(pred,actual)={corr:.6f}, "
              f"rel_rmse={rel_rmse:.4f}, E[min_std_N(0,1)]={c_n3:.6f}")

print()
print("=> High correlation = Gaussian approximation is good.")
print("=> If Gaussian, then m2m = 1 - |c_n3| * CoV, and higher CoV => lower m2m.")
print("=> This would close the proof: c2_c0 < R follows from CoV(v2-col) > CoV(v0-col).")
print()
print("done")
