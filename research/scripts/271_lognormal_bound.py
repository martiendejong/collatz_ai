"""
271_lognormal_bound.py
======================
FORMAL PROOF ATTEMPT for step (3b) via log-normal approximation bound.

STRATEGY:
  For a log-normal triplet (X1,X2,X3) with iid log Xi ~ N(mu, sig^2):
    m2m_LN(sig) = E[min(X)] / E[mean(X)] = 3*E[e^(mu+sig*Z1) * 1_{Z1<=Z2,Z1<=Z3}] / e^(mu+sig^2/2)
  This is ANALYTICALLY COMPUTABLE and strictly decreasing in sig.

  For K-L column triplets: if we can show
    (i) The K-L column log-distribution is "close to" log-normal with sigma=sigma_v0,sigma_v2
    (ii) sigma_v2 > sigma_v0 (follows from CoV^2(v2) > CoV^2(v0))
    (iii) The log-normal m2m function is monotone decreasing with margin epsilon
    (iv) The approximation error |m2m_KL - m2m_LN| < epsilon/2
  Then: m2m_v2 < m2m_v0 (proved for all k).

  KEY FACT: For iid log-normal triplets (X1,X2,X3) with log Xi ~ N(0, sig^2):
    m2m_LN(sig) = E[min(e^(sig*Z1), e^(sig*Z2), e^(sig*Z3))] / E[e^(sig*Z)]
               = E[min(e^(sig*Zi))] / e^(sig^2/2)
    Using order statistics: E[X_(1)] for n=3 iid log-normals.

  For ordered std normals Z_(1) <= Z_(2) <= Z_(3):
    E[Z_(1)] = -3/sqrt(pi) * int_0^inf Phi(-z)^2 * phi(z) dz = -3/sqrt(2*pi) * E[|Z|*Phi(-Z)^2]
  For n=3 iid N(0,1): E[Z_(1)] = -3/sqrt(pi) * Gamma(1+1/2)/(1+1/2)... [standard result]:
    E[Z_(1)] ≈ -0.8463 for n=3.

  m2m_LN(sig) = exp(sig * E[Z_(1)]) = exp(-0.8463 * sig)   [for iid symmetric case]

  This is STRICTLY DECREASING in sig > 0.

PLAN:
  1. Compute log-sigma for v0 and v2 column triplets.
  2. Verify sigma_v2 > sigma_v0 (consistent with CoV^2 ordering).
  3. Compute m2m_LN for each column and compare to actual m2m.
  4. Bound the approximation error.
  5. Check: gap(m2m_v0 - m2m_v2) > 2 * max_error => proof holds.
"""
import numpy as np
from math import log2, sqrt, exp, pi
from scipy import stats, integrate
import sys

ALPHA = log2(3.0)

def run_kl(k, lam, n_iter=None):
    if n_iter is None:
        n_iter = 800 + 200*max(0, k-8)
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

# ======================================================================
# Log-normal m2m function: E[min of 3 iid LN] / E[LN]
# For iid N(0, sig^2): E[min] = e^(sig^2/2) * E[min of 3 iid LN(0,sig^2)]
# = e^(sig^2/2) * e^(sig * E[Z_(1)])  where E[Z_(1)] is expected min of 3 std normals
# Actually: if X_i = exp(sig * Z_i), Z_i iid N(0,1):
#   E[X_(1)] = E[min(exp(sig*Z_i))] = E[exp(sig*Z_(1))]
#   E[X] = E[exp(sig*Z)] = exp(sig^2/2)
#   m2m_LN = E[exp(sig*Z_(1))] / exp(sig^2/2)
#
# For n=3 iid Z~N(0,1):
#   f_{Z_(1)}(z) = 3 * phi(z) * Phi(z)^2   [Wait: Z_(1) = min, so F_{Z_(1)}(z) = 1-(1-Phi(z))^3]
#   f_{Z_(1)}(z) = 3 * phi(z) * (1-Phi(z))^2
# E[exp(sig*Z_(1))] = 3 * int_{-inf}^{inf} exp(sig*z) * phi(z) * (1-Phi(z))^2 dz
#
# Numerical computation:

def m2m_lognormal(sig, n_samples=1000000):
    """E[min of 3 iid LN(0,sig^2)] / E[LN(0,sig^2)] via Monte Carlo."""
    if sig < 1e-10:
        return 1.0
    rng = np.random.default_rng(42)
    Z = rng.normal(0, sig, size=(n_samples, 3))
    X = np.exp(Z)
    m2m = X.min(axis=1) / X.mean(axis=1)
    return float(np.mean(m2m))

def m2m_lognormal_analytic(sig):
    """E[min of 3 iid LN(0,sig^2)] / E[LN(0,sig^2)] analytically.
    = 3 * int exp(sig*z) * phi(z) * (1-Phi(z))^2 dz / exp(sig^2/2)
    """
    if sig < 1e-10:
        return 1.0
    def integrand(z):
        return np.exp(sig*z) * stats.norm.pdf(z) * (1 - stats.norm.cdf(z))**2
    I, _ = integrate.quad(integrand, -20, 20, limit=200)
    return 3 * I / np.exp(sig**2 / 2)

# Verify: m2m_LN is strictly decreasing in sig
print("Log-normal m2m function (strictly decreasing in sig):")
print(f"{'sig':>8} {'m2m_MC':>10} {'m2m_analytic':>14}")
sigs = [0.01, 0.02, 0.03, 0.05, 0.08, 0.10, 0.15, 0.20]
for sig in sigs:
    mc = m2m_lognormal(sig, 500000)
    an = m2m_lognormal_analytic(sig)
    print(f"sig={sig:.3f}: m2m_MC={mc:.7f}, m2m_analytic={an:.7f}")
sys.stdout.flush()

# ======================================================================
# Main analysis: K-L column triplets vs log-normal
# ======================================================================
def analyze(k, lam):
    v, Nl, A, B1, B3 = run_kl(k, lam)
    Nl3 = Nl // 3
    j3  = np.arange(Nl3, dtype=np.int64)
    v0 = v[0::3]; v2 = v[2::3]

    # v0 and v2 column triplets
    col0 = np.stack([v0[j3], v0[j3+Nl3], v0[j3+2*Nl3]], axis=1)  # (Nl3, 3)
    col2 = np.stack([v2[j3], v2[j3+Nl3], v2[j3+2*Nl3]], axis=1)

    # Log-columns
    logcol0 = np.log(col0 + 1e-300)
    logcol2 = np.log(col2 + 1e-300)

    # Within-column log-std (proxy for sigma of log-normal model)
    logsig0 = logcol0.std(axis=1)   # within-column log-std per triplet
    logsig2 = logcol2.std(axis=1)

    mean_logsig0 = float(np.mean(logsig0))
    mean_logsig2 = float(np.mean(logsig2))

    # Per-column m2m ratio
    mean_col0 = col0.mean(axis=1)
    mean_col2 = col2.mean(axis=1)
    min_col0  = col0.min(axis=1)
    min_col2  = col2.min(axis=1)
    mmr_col0  = min_col0 / (mean_col0 + 1e-30)
    mmr_col2  = min_col2 / (mean_col2 + 1e-30)

    # Log-normal prediction for each column
    # m2m_LN(sig) ≈ exp(-0.8463 * sig) for iid symmetric case
    # (0.8463 = |E[min of 3 iid N(0,1)]|)
    # Use the analytic formula for each column's sigma
    EZ1 = -0.84628  # E[min of 3 iid N(0,1)], computed analytically

    m2m_LN_pred0 = np.exp(EZ1 * logsig0)
    m2m_LN_pred2 = np.exp(EZ1 * logsig2)

    # Approximation error per column
    err0 = np.abs(mmr_col0 - m2m_LN_pred0)
    err2 = np.abs(mmr_col2 - m2m_LN_pred2)

    max_err0 = float(np.max(err0))
    max_err2 = float(np.max(err2))
    mean_err0 = float(np.mean(err0))
    mean_err2 = float(np.mean(err2))

    # Overall m2m
    m2m_v0 = float(np.mean(mmr_col0))
    m2m_v2 = float(np.mean(mmr_col2))
    gap_m2m = m2m_v0 - m2m_v2

    # Log-normal prediction for overall m2m
    m2m_LN_v0 = float(np.mean(m2m_LN_pred0))
    m2m_LN_v2 = float(np.mean(m2m_LN_pred2))
    gap_LN = m2m_LN_v0 - m2m_LN_v2

    # Is gap > 2 * max_error?
    total_max_err = max_err0 + max_err2
    proof_ok = gap_m2m > total_max_err

    # rho and R
    s = np.arange(Nl, dtype=np.int64)
    sigma1 = (4*s+2) % Nl
    v1 = v[1::3]
    rho = A / float(np.mean(v1 / v0[sigma1]))
    t   = A / rho
    R_th = (t**2 + lam) / (1 + t*lam)
    mean_v0 = float(np.mean(v0))
    mean_v2 = float(np.mean(v2))
    c0 = float(np.mean(min_col0))
    c2 = float(np.mean(min_col2))

    return {
        'k': k, 'lam': lam, 't': t, 'R': R_th,
        'c0': c0, 'c2': c2, 'c2/c0': c2/c0, 'c2/c0<=R': c2/c0 <= R_th,
        'm2m_v0': m2m_v0, 'm2m_v2': m2m_v2, 'gap_m2m': gap_m2m,
        'mean_logsig0': mean_logsig0, 'mean_logsig2': mean_logsig2,
        'm2m_LN_v0': m2m_LN_v0, 'm2m_LN_v2': m2m_LN_v2, 'gap_LN': gap_LN,
        'max_err0': max_err0, 'max_err2': max_err2, 'total_max_err': total_max_err,
        'mean_err0': mean_err0, 'mean_err2': mean_err2,
        'proof_ok': proof_ok,
        'gap_vs_err': gap_m2m / total_max_err,
    }

print("\n" + "="*70)
print("K-L COLUMN TRIPLET: LOG-NORMAL APPROXIMATION QUALITY")
print(f"{'lam':>6} {'k':>3} {'sig0':>7} {'sig2':>7} {'m2m_v0':>8} {'m2m_v2':>8} {'gap':>7} {'maxerr':>8} {'ratio':>7} {'ok':>5}")
for lam in [1.30, 1.50, 1.70, 1.90]:
    for k in [8, 10]:
        d = analyze(k, lam)
        print(f"lam={lam:.2f} k={k:>2} {d['mean_logsig0']:>7.5f} {d['mean_logsig2']:>7.5f} {d['m2m_v0']:>8.5f} {d['m2m_v2']:>8.5f} {d['gap_m2m']:>7.5f} {d['total_max_err']:>8.5f} {d['gap_vs_err']:>7.2f} {str(d['proof_ok']):>5}")
    sys.stdout.flush()

print()
# Detailed at lam=1.70, k=10
lam, k = 1.70, 10
d = analyze(k, lam)
print(f"\nDetailed: lam={lam}, k={k}")
print(f"  Mean log-sigma v0: {d['mean_logsig0']:.6f}, v2: {d['mean_logsig2']:.6f}")
print(f"  LN prediction m2m v0: {d['m2m_LN_v0']:.6f}, v2: {d['m2m_LN_v2']:.6f}")
print(f"  Actual m2m v0: {d['m2m_v0']:.6f}, v2: {d['m2m_v2']:.6f}")
print(f"  Gap (m2m_v0 - m2m_v2): {d['gap_m2m']:.6f}")
print(f"  Max LN approx error v0: {d['max_err0']:.6f}, v2: {d['max_err2']:.6f}")
print(f"  Total max error: {d['total_max_err']:.6f}")
print(f"  Gap / total_error: {d['gap_vs_err']:.2f}")
print(f"  Proof via LN bound: {d['proof_ok']}")
print()

print("""
PROOF STRUCTURE (step 3b via log-normal bound):

For each K-L column triplet (X0, X1, X2) at level k:
  1. Fit log-sigma: sig_j = std(log(X0_j, X1_j, X2_j))
  2. Log-normal prediction: m2m_LN(j) = exp(EZ1 * sig_j) where EZ1 = -0.84628
  3. Actual: m2m(j) = min(X_j) / mean(X_j)
  4. Error: |m2m(j) - m2m_LN(j)| <= epsilon_max

If epsilon_max for v0 and v2 together is smaller than gap(m2m_v0, m2m_v2):
  => gap(m2m_v0, m2m_v2) > 0 is PROVED.

Since m2m_LN is strictly decreasing and mean_logsig_v2 > mean_logsig_v0:
  m2m_LN_v2 < m2m_LN_v0 (log-normal bound).
If |m2m - m2m_LN| <= epsilon_max for all columns:
  m2m_v2 <= m2m_LN_v2 + epsilon_max < m2m_LN_v0 - (gap_LN) + epsilon_max
  For this to give m2m_v2 < m2m_v0: need gap_LN > 2*epsilon_max.

FINDING: Whether gap > 2*max_error determines if the LN-bound proof closes step (3b).
""")
print("done")
