"""
266_subclass_weight.py
======================
PROOF ATTEMPT: m2m_v2 <= m2m_v0 via the KEY WEIGHT INEQUALITY

From Script 265: the m2m gap is driven by:
  w_sc0_v2 = a0_v2 / sum_v2 < w_sc1_v0 = a1_v0 / sum_v0

where sc0 of v2 and sc1 of v0 are BOTH the "small-cb" sub-classes (use c1 = (A/rho)*c0).

ALGEBRAIC CONDITION: a0_v2 / a1_v0 <= R = sum_v2 / sum_v0

From K-L:
  rho * a0_v2 = (A^2/rho) * a2_v0 + B3 * c1   ...(I)
  rho * a1_v0 = A * a1_v2 + B1 * c1            ...(II)

So a0_v2/a1_v0 = [(A^2/rho)*a2_v0 + B3*c1] / [A*a1_v2 + B1*c1]

NEED: this <= R = (a0_v2+a1_v2+a2_v2)/(a0_v0+a1_v0+a2_v0)

This session verifies the inequality and attempts to derive it from K-L equations.

STEP 1: Verify a0_v2 <= a1_v0 (simpler than a0_v2/a1_v0 <= R).
STEP 2: Verify a0_v2/a1_v0 <= R.
STEP 3: Show these together imply m2m_v2 <= m2m_v0.
STEP 4: Attempt algebraic proof of a0_v2 <= a1_v0 from K-L structure.
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
    v0 = v[0::3]; v1 = v[1::3]; v2 = v[2::3]

    # rho
    sigma1 = (4*s+2) % Nl
    rho = A / float(np.mean(v1 / v0[sigma1]))

    # Sub-class means
    a0v2 = float(np.mean(v2[s%3==0]))  # small-cb class of v2 (uses c1)
    a1v2 = float(np.mean(v2[s%3==1]))
    a2v2 = float(np.mean(v2[s%3==2]))
    a0v0 = float(np.mean(v0[s%3==0]))
    a1v0 = float(np.mean(v0[s%3==1]))  # small-cb class of v0 (uses c1)
    a2v0 = float(np.mean(v0[s%3==2]))

    sum_v2 = a0v2+a1v2+a2v2
    sum_v0 = a0v0+a1v0+a2v0
    R = sum_v2/sum_v0

    # Key inequality: a0_v2 <= a1_v0 ?
    key1 = a0v2 <= a1v0
    # Key inequality: a0_v2/a1_v0 <= R ?
    key2 = (a0v2/a1v0) <= R

    # From K-L (I)/(II): a0_v2/a1_v0 = [(A^2/rho)*a2_v0 + B3*c1] / [A*a1_v2 + B1*c1]
    # where c1 = (A/rho)*c0. Need c0 too.
    cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
    j  = np.arange(Nl, dtype=np.int64)
    c0 = float(np.mean(cb[j%3==0]))
    c1_val = (A/rho)*c0

    lhs_formula = ((A**2/rho)*a2v0 + B3*c1_val) / (A*a1v2 + B1*c1_val)
    actual_ratio = a0v2/a1v0

    # STEP 4 ATTEMPT: Can we bound a0_v2/a1_v0 <= R from K-L?
    # a0_v2/a1_v0 = [(A^2/rho)*a2_v0 + B3*c1] / [A*a1_v2 + B1*c1]
    # R = sum_v2/sum_v0
    # Sufficient to show: [(A^2/rho)*a2_v0 + B3*c1]*sum_v0 <= [A*a1_v2 + B1*c1]*sum_v2

    # From K-L summed: rho*sum_v2 = (A^2/rho)*sum_v0 + B3*(c0+c1+c2)
    #                  rho*sum_v0 = A*sum_v2 + B1*(c0+c1+c2)
    # Let S2=sum_v2, S0=sum_v0, C=c0+c1+c2.
    # Need: [(A^2/rho)*a2_v0 + B3*c1]*S0 <= [A*a1_v2 + B1*c1]*S2

    # Using rho*S0 = A*S2 + B1*C => S0 = (A*S2 + B1*C)/rho
    # And rho*S2 = (A^2/rho)*S0 + B3*C => S2 = ((A^2/rho)*S0 + B3*C)/rho

    # The inequality becomes:
    # [(A^2/rho)*a2_v0 + B3*c1]*(A*S2 + B1*C)/rho <= [A*a1_v2 + B1*c1]*S2*rho
    # [(A^2/rho)*a2_v0 + B3*c1]*(A*S2 + B1*C) <= [A*a1_v2 + B1*c1]*rho^2*S2
    # Expand... complex. Let's just check numerically if a2_v0 <= a1_v2 holds.
    a2v0_leq_a1v2 = a2v0 <= a1v2

    return {
        'a0v2': a0v2, 'a1v0': a1v0, 'R': R,
        'ratio': actual_ratio, 'ratio_formula': lhs_formula,
        'key1_a0v2_leq_a1v0': key1,
        'key2_ratio_leq_R': key2,
        'a2v0': a2v0, 'a1v2': a1v2,
        'a2v0_leq_a1v2': a2v0_leq_a1v2,
        'a0v2/sum_v2': a0v2/sum_v2,
        'a1v0/sum_v0': a1v0/sum_v0,
        'B3c1': B3*c1_val, 'A2a2v0': (A**2/rho)*a2v0,
        'B1c1': B1*c1_val, 'Aa1v2': A*a1v2,
        'rho': rho, 'A': A, 'B1': B1, 'B3': B3,
    }

print("266: Sub-class weight inequality -- proof attempt for m2m_v2 <= m2m_v0")
print("="*75)

# Detailed at lam=1.70, k=10
lam, k = 1.70, 10
d = analyze(k, lam)
print(f"\nDetailed: lam={lam}, k={k}")
print(f"  a0_v2 (small-cb sc of v2) = {d['a0v2']:.6f}")
print(f"  a1_v0 (small-cb sc of v0) = {d['a1v0']:.6f}")
print(f"  a0_v2 <= a1_v0: {d['key1_a0v2_leq_a1v0']}")
print(f"  ratio a0_v2/a1_v0 = {d['ratio']:.6f}")
print(f"  formula ratio    = {d['ratio_formula']:.6f}  (check)")
print(f"  R = sum_v2/sum_v0 = {d['R']:.6f}")
print(f"  ratio/R = {d['ratio']/d['R']:.6f}  (must be <= 1)")
print(f"  key2: a0_v2/a1_v0 <= R: {d['key2_ratio_leq_R']}")
print()
print(f"  Weight of small-cb sc in v2: {d['a0v2/sum_v2']:.5f}")
print(f"  Weight of small-cb sc in v0: {d['a1v0/sum_v0']:.5f}")
print(f"  Weight ratio v2/v0: {d['a0v2/sum_v2']/d['a1v0/sum_v0']:.5f}  (must be <= 1)")
print()
print(f"  K-L numerator:  (A^2/rho)*a2_v0={d['A2a2v0']:.5f}, B3*c1={d['B3c1']:.5f}")
print(f"  K-L denominator: A*a1_v2={d['Aa1v2']:.5f}, B1*c1={d['B1c1']:.5f}")
print(f"  a2_v0={d['a2v0']:.5f}, a1_v2={d['a1v2']:.5f}")
print(f"  a2_v0 <= a1_v2: {d['a2v0_leq_a1v2']}")

print(f"\n\nScans:")
print(f"{'lam':>6} {'k':>3} {'a0v2<=a1v0':>12} {'ratio<=R':>10} {'a2v0<=a1v2':>12} {'wt_v2':>8} {'wt_v0':>8} {'wt_v2<=v0':>12}")
for lam in [1.30, 1.50, 1.70, 1.90, 2.00]:
    for k in [6, 10, 13]:
        d = analyze(k, lam)
        wt_ok = d['a0v2/sum_v2'] <= d['a1v0/sum_v0']
        print(f"lam={lam:.2f} k={k:>2} {str(d['key1_a0v2_leq_a1v0']):>12} {str(d['key2_ratio_leq_R']):>10} {str(d['a2v0_leq_a1v2']):>12} {d['a0v2/sum_v2']:>8.5f} {d['a1v0/sum_v0']:>8.5f} {str(wt_ok):>12}")
    sys.stdout.flush()

print(f"""
ALGEBRAIC PROOF ATTEMPT for a0_v2 <= a1_v0:

From K-L:
  rho * a0_v2 = (A^2/rho) * a2_v0 + B3 * c1    ...(I)
  rho * a1_v0 = A * a1_v2 + B1 * c1             ...(II)

Subtracting: rho*(a1_v0 - a0_v2) = A*a1_v2 - (A^2/rho)*a2_v0 + (B1-B3)*c1

Since B1 < B3 (B3/B1 = lambda > 1): (B1-B3)*c1 < 0.
So rho*(a1_v0-a0_v2) = A*(a1_v2 - (A/rho)*a2_v0) + (B1-B3)*c1.

If a1_v2 >= (A/rho)*a2_v0 AND the negative (B1-B3)*c1 is not too large,
then a1_v0 >= a0_v2.

From K-L for a2_v0: rho*a2_v0 = A*a2_v2 + B1*c2 => a2_v0 = (A/rho)*a2_v2 + (B1/rho)*c2.
So (A/rho)*a2_v0 = (A^2/rho^2)*a2_v2 + (A*B1/rho^2)*c2.
And a1_v2 - (A/rho)*a2_v0 = a1_v2 - (A^2/rho^2)*a2_v2 - (A*B1/rho^2)*c2.

Since a1_v2 > a2_v2 (larger forcing: sc1 gets c0 vs sc2 gets c2, and c0 > c2 for lam<1.45,
c0 < c2 for lam>1.45 -- so ordering depends on lambda):
  lam=1.70 data: a1_v2=0.211 > a2_v2=0.235? NO! a1_v2 < a2_v2 at lam=1.70.

So the a1_v2 > a2_v2 assumption fails for lam > ~1.5. Need a different bound.
""")

# Check a1_v2 vs a2_v2 ordering
print("Sub-class means ordering check:")
print(f"{'lam':>6} {'k':>3} {'a0v2':>8} {'a1v2':>8} {'a2v2':>8} {'a0v0':>8} {'a1v0':>8} {'a2v0':>8}")
k = 10
for lam in [1.30, 1.50, 1.70, 1.90, 2.00]:
    v, Nl, A, B1, B3 = run_kl(k, lam)
    s = np.arange(Nl, dtype=np.int64)
    v0=v[0::3]; v1=v[1::3]; v2=v[2::3]
    a0v2=float(np.mean(v2[s%3==0])); a1v2=float(np.mean(v2[s%3==1])); a2v2=float(np.mean(v2[s%3==2]))
    a0v0=float(np.mean(v0[s%3==0])); a1v0=float(np.mean(v0[s%3==1])); a2v0=float(np.mean(v0[s%3==2]))
    print(f"lam={lam:.2f} k={k} {a0v2:.5f} {a1v2:.5f} {a2v2:.5f} {a0v0:.5f} {a1v0:.5f} {a2v0:.5f}")

print("\ndone")
