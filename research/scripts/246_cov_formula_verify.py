"""
246_cov_formula_verify.py
=========================
Verify the analytical formula Cov(X, min(X,Y,Z)) = (1/3)(E[min^2] - mu*E[min])
for the K-L eigenvector CODE-triplets.

Also verify the negativity condition: (mu - E[min])*E[min] > Var[min]
at all depths k and lambda values.

The formula:
  Cov(X, min(X,Y,Z)) = E[X*min] - E[X]*E[min]
  For exchangeable (X,Y,Z): 3*Cov = Cov(X+Y+Z, min) = E[(X+Y+Z)*min] - 3*mu*E[min]
  = E[sum*min] - 3*mu*E[min]

  Exact derivation for exchangeable:
    E[X*min(X,Y,Z)] = E[X*min | X=min]*P(X=min) + E[X*min(Y,Z) | X>min(Y,Z)]*P(X>min(Y,Z))
    = (1/3)*E[X^2 | X=min] + (2/3)*E[X]*E[min(Y,Z)]
    = (1/3)*E[min^2] + (2/3)*mu*E[min]
    (using X and min(Y,Z) independent when X>min(Y,Z), E[X]=mu, E[min(Y,Z)]=E[min])

  Therefore: Cov(X, min(X,Y,Z)) = (1/3)*E[min^2] + (2/3)*mu*E[min] - mu*E[min]
    = (1/3)*(E[min^2] - mu*E[min])

Note: The formula uses approximate independence of X from min(Y,Z) when X>min(Y,Z),
which is exact for i.i.d. and approximate for exchangeable correlated triplets.
"""
import sys
import numpy as np
from math import log2

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

def analyze_triplets(v, Nl):
    """Analyze the CODE-triplets of v2 and check the Cov formula."""
    Nl3 = Nl // 3
    s0 = np.arange(Nl3, dtype=np.int64)
    v2 = v[2::3]  # r=2 component, length Nl

    # For each j-equiv-2 mod 3 group at the CODE-triplet level, we need the
    # sigma0 pullback. But simpler: look at the v2 CODE-triplets directly.
    # The j=2 group (for k=4, Nl=9, Nl3=3) uses v2 values {8, 2, 5}.
    # More generally: the CODE-triplet {j, j+Nl3, j+2Nl3} with sigma0 pullback.

    # The sigma0 map permutes v2 values in CODE-triplets.
    # Let's directly analyze the v2 CODE-triplets (before sigma0 pullback).
    # Each CODE-triplet = (v2[s], v2[s+Nl3], v2[s+2Nl3]) for s=0..Nl3-1.

    a = v2[s0]          # shape (Nl3,)
    b = v2[s0 + Nl3]
    c = v2[s0 + 2*Nl3]
    min_abc = np.minimum(np.minimum(a, b), c)

    mu = float(np.mean(a))  # = mean of v2 (all Nl values, approx)
    E_min = float(np.mean(min_abc))
    E_min2 = float(np.mean(min_abc**2))
    Var_min = float(np.var(min_abc))

    # Formula: Cov(a, min(a,b,c)) = (1/3)*(E[min^2] - mu*E[min])
    cov_formula = (E_min2 - mu * E_min) / 3.0

    # Direct measurement: Cov(a, min(a,b,c))
    cov_direct = float(np.mean(a * min_abc)) - mu * E_min

    # Negativity condition: (mu - E[min])*E[min] > Var[min]
    lhs = (mu - E_min) * E_min
    rhs = Var_min
    neg_cond = lhs > rhs

    return {
        'mu': mu, 'E_min': E_min, 'E_min2': E_min2, 'Var_min': Var_min,
        'cov_formula': cov_formula, 'cov_direct': cov_direct, 'neg_cond': neg_cond,
        'lhs': lhs, 'rhs': rhs
    }

print("246: Verify Cov(X, min(X,Y,Z)) = (1/3)(E[min^2] - mu*E[min])")
print("Also verify negativity condition: (mu-E[min])*E[min] > Var[min]")
print("="*80)

# Lambda scan at k=12
K = 12
LAMS = [1.30, 1.50, 1.70, 1.90, 2.00]
print(f"\nLambda scan at k={K}:")
print(f"{'lam':>6}  {'cov_formula':>13}  {'cov_direct':>13}  {'err':>8}  {'(mu-Em)*Em':>12}  {'Var_min':>10}  {'neg?':>5}")
for lam in LAMS:
    v, Nl = run_kl(K, lam)
    res = analyze_triplets(v, Nl)
    err = abs(res['cov_formula'] - res['cov_direct']) / max(abs(res['cov_direct']), 1e-15)
    print(f"lam={lam:.2f}  {res['cov_formula']:>13.6e}  {res['cov_direct']:>13.6e}  "
          f"{err:>8.2e}  {res['lhs']:>12.6e}  {res['rhs']:>10.6e}  {'YES' if res['neg_cond'] else 'NO':>5}")
    sys.stdout.flush()

print()
# Depth scan at lambda=1.70
LAM = 1.70
print(f"Depth scan at lambda={LAM}:")
print(f"{'k':>4}  {'cov_formula':>13}  {'cov_direct':>13}  {'err':>8}  {'(mu-Em)*Em':>12}  {'Var_min':>10}  {'neg?':>5}")
for k in range(4, 14):
    v, Nl = run_kl(k, LAM)
    res = analyze_triplets(v, Nl)
    err = abs(res['cov_formula'] - res['cov_direct']) / max(abs(res['cov_direct']), 1e-15)
    print(f"k={k:>2}  {res['cov_formula']:>13.6e}  {res['cov_direct']:>13.6e}  "
          f"{err:>8.2e}  {res['lhs']:>12.6e}  {res['rhs']:>10.6e}  {'YES' if res['neg_cond'] else 'NO':>5}")
    sys.stdout.flush()

print()
print("Note: Formula is exact for i.i.d., approximate for correlated triplets.")
print("Error measures departure from i.i.d. assumption due to CODE-triplet correlations.")
print()
print("done")
