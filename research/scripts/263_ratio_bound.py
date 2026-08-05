"""
263_ratio_bound.py
==================
ANALYTICAL CLOSURE ATTEMPT: Can we bound c2/c0 via mean_v2/mean_v0?

From the K-L scalar mean equations:
  rho * mean_v2 = (A^2/rho) * mean_v0 + B3 * cbar
  rho * mean_v0 = A * mean_v2 + B1 * cbar

Solving: R = mean_v2/mean_v0 = (A^2/rho^2 + lambda) / (1 + A*lambda/rho)
                               = (A^2 + rho^2*lambda) / (rho^2 + A*rho*lambda)

KEY QUESTIONS:
  Q1: Does c2/c0 <= mean_v2/mean_v0 = R(lambda) hold numerically?
  Q2: If yes, is R(lambda) << 1 + RHS(lambda) ?
  Q3: Can R(lambda) - 1 < RHS(lambda) be proved analytically?
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
    cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
    j  = np.arange(Nl, dtype=np.int64)
    c0 = float(np.mean(cb[j%3==0]))
    c2 = float(np.mean(cb[j%3==2]))

    # Eigenvector type means
    r_arr = np.arange(3**(k-1), dtype=np.int64) % 3
    mean_v0 = float(np.mean(v[r_arr==0]))  # r=0 type nodes
    mean_v2 = float(np.mean(v[r_arr==2]))  # r=2 type nodes

    # rho estimate (from r=1 equation)
    v1 = v[1::3]; v0 = v[0::3]
    s  = np.arange(Nl, dtype=np.int64)
    sigma1 = (4*s+2) % Nl
    rho = A / float(np.mean(v1 / v0[sigma1]))

    # Analytical R formula
    R_analytical = (A**2/rho**2 + lam) / (1 + A*lam/rho)

    # RHS of the gap
    rhs = lam**5 * (1 - A/rho) * rho**2

    return c0, c2, mean_v0, mean_v2, rho, A, R_analytical, rhs

print("263: Bound c2/c0 via mean_v2/mean_v0 = R(lambda)")
print("="*80)

# Q1: Is c2/c0 <= mean_v2/mean_v0 ?
print(f"\nScan lam=1.70, k=5..13:")
print(f"{'k':>4}  {'c2/c0':>8}  {'R_num':>8}  {'R_anal':>8}  {'c2/c0<=R?':>10}  {'(R-1)':>8}  {'RHS':>8}  {'mar':>6}")
lam = 1.70
for k in range(5, 14):
    c0, c2, mv0, mv2, rho, A, R_anal, rhs = analyze(k, lam)
    ratio = c2/c0
    R_num  = mv2/mv0
    ok = ratio <= R_num
    print(f"k={k:>2}  {ratio:>8.5f}  {R_num:>8.5f}  {R_anal:>8.5f}  {str(ok):>10}  {R_num-1:>8.5f}  {rhs:>8.4f}  {rhs/(R_num-1):>6.1f}")
    sys.stdout.flush()

print(f"\nLambda scan k=10:")
print(f"{'lam':>6}  {'c2/c0':>8}  {'R_num':>8}  {'R_anal':>8}  {'c2/c0<=R?':>10}  {'(R-1)':>8}  {'RHS':>8}  {'mar':>6}")
k = 10
for lam in [1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00]:
    c0, c2, mv0, mv2, rho, A, R_anal, rhs = analyze(k, lam)
    ratio = c2/c0
    R_num  = mv2/mv0
    ok = ratio <= R_num
    excess = max(R_num - 1, 1e-8)
    print(f"lam={lam:.2f}  {ratio:>8.5f}  {R_num:>8.5f}  {R_anal:>8.5f}  {str(ok):>10}  {R_num-1:>8.5f}  {rhs:>8.4f}  {rhs/excess:>6.1f}")
    sys.stdout.flush()

# Q3: Prove R(lambda) - 1 < RHS analytically
print(f"""
Q3 ANALYTICAL: Is R(lambda)-1 < RHS(lambda) = lambda^5*(1-A/rho)*rho^2 ?

R(lam) = (A^2/rho^2 + lam) / (1 + A*lam/rho)
R - 1  = (A^2/rho^2 + lam - 1 - A*lam/rho) / (1 + A*lam/rho)
       = ((A/rho)^2 - (A/rho)*lam + lam - 1) / (1 + A*lam/rho)

Let t = A/rho (in (0,1) for lam>1).
R - 1 = (t^2 - t*lam + lam - 1) / (1 + t*lam)
      = (lam*(1-t) - (1-t^2)) / (1 + t*lam)
      = (1-t)*(lam - 1 - t) / (1 + t*lam)   [IMPORTANT: = (1-t)*(lam-1-t)/(1+t*lam)]

RHS = lam^5 * (1-t) * rho^2

So: (R-1)/RHS = (lam-1-t) / (lam^5 * rho^2 * (1 + t*lam))

For lam>1, t<1: lam-1-t = lam-1-A/rho. Since A=lam^-2, rho>1 for lam>1:
  t = A/rho = lam^-2/rho.
  lam-1-t = lam-1-lam^-2/rho ~ lam-1 for large lam (t->0).

Denominator: lam^5 * rho^2 * (1+t*lam) ~ lam^5 for large lam.

So (R-1)/RHS ~ (lam-1) / lam^5 -> 0 as lam -> inf.

For lam near 1: lam-1 ~ eps, t ~ 1-O(eps), so lam-1-t ~ eps-(1-O(eps)) < 0!
  => R-1 <= 0 for lam near 1 (c2 <= c0), consistent with data (lambda<1.45).

CONCLUSION OF Q3:
  For lam < some lam_0 (approx 1.45): R-1 < 0 => R < 1 => c2/c0 <= R < 1 <= 1+RHS. OK.
  For lam > lam_0: (R-1)/RHS = (lam-1-t)/(lam^5*rho^2*(1+t*lam)) is small and
    decreasing in lam (numerator linear, denominator lam^5). Maximum around lam~1.5-1.8.

So R(lam) - 1 < RHS(lam) holds for ALL lam > 1 if c2/c0 <= R(lam).
""")
sys.stdout.flush()

# Verify the factored form (R-1)/RHS numerically
print("Numerical (R-1)/RHS ratio (must be < 1 for all lam,k):")
print(f"{'lam':>6}  {'k':>3}  {'R-1':>8}  {'RHS':>8}  {'(R-1)/RHS':>10}  {'c2/c0<=R':>10}")
for lam in [1.30, 1.50, 1.70, 1.90, 2.00]:
    k = 10
    c0, c2, mv0, mv2, rho, A, R_anal, rhs = analyze(k, lam)
    R_num = mv2/mv0
    t = A/rho
    r_minus_1 = R_num - 1
    rat = r_minus_1 / rhs if abs(rhs) > 1e-10 else float('inf')
    ok = c2/c0 <= R_num
    print(f"lam={lam:.2f}  k={k}  {r_minus_1:>8.4f}  {rhs:>8.4f}  {rat:>10.4f}  {str(ok):>10}")

print("\ndone")
