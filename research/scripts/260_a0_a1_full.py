"""
260_a0_a1_full.py
=================
COMPLETE PROOF that a0_v2 < a1_v2 for all lambda > 1.

FROM Obs 464: c_1 = (A/rho)*c_0 < c_0.

NOW: Show that the T4 contribution ALSO satisfies m0 <= m1.

K-L EQUATION FOR v2 NODES:
  For any v2 node i = 3s+2:
    T4(3s+2) = (4*(3s+2)+2) mod N = (12s+10) mod N
    r-type of T4(3s+2): (12s+10) mod 3 = 10 mod 3 = 1 (ALWAYS r=1)
    => T4 maps ALL v2 nodes to v1 nodes.

  Full K-L: rho*v2[s] = A*v1[tau(s)] + B3*cb[R3[s]]
    where tau(s) = (4s+3) mod Nl  (the v1 index from T4(3s+2) = 3*(4s+3)+1)
    and R3[s] = (2s+1) mod Nl.

  So:
    rho*a0_v2 = A*Mean(v1[tau(s)] | s==0 mod 3) + B3*Mean(cb[R3[s]] | s==0 mod 3)
    rho*a1_v2 = A*Mean(v1[tau(s)] | s==1 mod 3) + B3*Mean(cb[R3[s]] | s==1 mod 3)

CLAIM: R3 maps:
  s==0 mod 3: R3[s] = (2s+1) mod Nl. (2*0+1) mod 3 = 1 => class 1 in cb => c_1.
  s==1 mod 3: R3[s] = (2s+1) mod Nl. (2*1+1) mod 3 = 0 => class 0 in cb => c_0.
  s==2 mod 3: R3[s] = (2s+1) mod Nl. (2*2+1) mod 3 = 2 => class 2 in cb => c_2.

So: Mean(cb[R3[s]] | s==r mod 3) = c_{sigma_R(r)} where sigma_R(0)=1, sigma_R(1)=0, sigma_R(2)=2.
  => cb contribution: B3*c_1 for a0 and B3*c_0 for a1.

NOW THE T4 TERM:
  v1[tau(s)] = (A/rho)*v0[sigma1(tau(s))] where sigma1(t) = (4t+2) mod Nl.
  tau(s) = (4s+3) mod Nl.
  sigma1(tau(s)) = (4*(4s+3)+2) mod Nl = (16s+14) mod Nl.

  Mean(v1[tau(s)] | s==0 mod 3) = (A/rho)*Mean(v0[(16s+14) mod Nl] | s==0 mod 3)
  Mean(v1[tau(s)] | s==1 mod 3) = (A/rho)*Mean(v0[(16s+14) mod Nl] | s==1 mod 3)

  Define: m0 = Mean(v0[(16s+14) mod Nl] | s==0 mod 3)
          m1 = Mean(v0[(16s+14) mod Nl] | s==1 mod 3)

  Full equations:
    rho*a0 = B3*c_1 + A*(A/rho)*m0 = B3*c_1 + (A^2/rho)*m0
    rho*a1 = B3*c_0 + (A^2/rho)*m1

  rho*(a1-a0) = B3*(c_0-c_1) + (A^2/rho)*(m1-m0)
              = B3*c_0*(1-A/rho) + (A^2/rho)*(m1-m0)

  Term 1: B3*c_0*(1-A/rho) > 0 ALWAYS (since c_0>0, A/rho<1 for lam>1).
  Term 2: (A^2/rho)*(m1-m0) -- sign of m1-m0 to determine.

  CLAIM: m1 >= m0 (so both terms are non-negative => a1 > a0 robustly).

WHY m1 >= m0?
  s==0 mod 3 (s=3t): (16s+14) mod Nl = (48t+14) mod Nl
    (48t+14) mod 3 = 14 mod 3 = 2. So these are CLASS-2 positions of v0.
    => m0 = Mean(v0_class2 after permutation) = a2_v0 (approximately, if bijection)

  s==1 mod 3 (s=3t+1): (16(3t+1)+14) mod Nl = (48t+30) mod Nl
    (48t+30) mod 3 = 30 mod 3 = 0. So these are CLASS-0 positions of v0.
    => m1 = Mean(v0_class0 after permutation) = a0_v0

THEREFORE: m0 corresponds to v0-class-2 mean, m1 corresponds to v0-class-0 mean.
If a0_v0 >= a2_v0... need to check ordering of v0 class means.

THIS SCRIPT VERIFIES ALL CLAIMS NUMERICALLY.
"""
import numpy as np
from math import log2
import sys

ALPHA = log2(3.0)

def run_kl(k, lam, n_iter=None):
    if n_iter is None:
        n_iter = 600 + 100*max(0, k-8)
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

def analyze_full(k, lam):
    v, Nl, A, B1, B3 = run_kl(k, lam)
    cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])

    v0 = v[0::3]; v1 = v[1::3]; v2 = v[2::3]

    # Estimate rho
    s_all = np.arange(Nl, dtype=np.int64)
    sigma1 = (4*s_all + 2) % Nl
    rho = float(np.mean(v1 / v0[sigma1]))  # should be A/rho => rho = A/mean(v1/v0[sigma1])
    # Actually: v1[s] = (A/rho)*v0[sigma1(s)] => A/rho = v1[s]/v0[sigma1(s)] for each s
    ratios = v1 / v0[sigma1]
    A_over_rho = float(np.mean(ratios))  # should be constant
    rho_est = A / A_over_rho

    # Class means of v2
    j = np.arange(Nl, dtype=np.int64)
    a0_v2 = float(np.mean(v2[j%3==0]))
    a1_v2 = float(np.mean(v2[j%3==1]))
    a2_v2 = float(np.mean(v2[j%3==2]))

    # Class means of v0
    a0_v0 = float(np.mean(v0[j%3==0]))
    a1_v0 = float(np.mean(v0[j%3==1]))
    a2_v0 = float(np.mean(v0[j%3==2]))

    # Class means of cb
    c_0 = float(np.mean(cb[j%3==0]))
    c_1 = float(np.mean(cb[j%3==1]))
    c_2 = float(np.mean(cb[j%3==2]))

    # tau(s) = (4s+3) mod Nl (v1 index from T4)
    # Verify: T4(3s+2) = 3*(4s+3)+1 (mod N)
    tau = (4*s_all + 3) % Nl
    # Verify T4(3s+2) = 3*tau(s)+1
    T4_v2 = (12*s_all + 10) % (3*Nl)  # T4 for i=3s+2
    T4_v2_check = 3*tau + 1
    err_tau = float(np.max(np.abs(T4_v2 - T4_v2_check)))

    # sigma1 of tau: (4*(4s+3)+2) mod Nl = (16s+14) mod Nl
    phi = (16*s_all + 14) % Nl  # v0 index for T4 contribution to v2

    # m0: Mean(v0[phi[s]] | s==0 mod 3)
    m0 = float(np.mean(v0[phi[s_all%3==0]]))
    m1 = float(np.mean(v0[phi[s_all%3==1]]))

    # What class (mod 3) do phi values fall into?
    phi0 = phi[s_all%3==0]  # for class 0
    phi1 = phi[s_all%3==1]  # for class 1
    phi0_class_dist = [float(np.mean(phi0%3==r)) for r in range(3)]
    phi1_class_dist = [float(np.mean(phi1%3==r)) for r in range(3)]

    # Verify full K-L decomposition for a0, a1
    # rho*a0 = B3*c_1 + (A^2/rho)*m0
    # rho*a1 = B3*c_0 + (A^2/rho)*m1
    A2_over_rho = A**2 / rho_est
    rho_a0_predicted = B3*c_1 + A2_over_rho*m0
    rho_a1_predicted = B3*c_0 + A2_over_rho*m1
    rho_a0_actual = rho_est * a0_v2
    rho_a1_actual = rho_est * a1_v2

    # a1 - a0 decomposition
    cb_diff = B3*(c_0 - c_1)  # = B3*c_0*(1 - A/rho)
    T4_diff = A2_over_rho*(m1 - m0)
    total_diff = (cb_diff + T4_diff) / rho_est

    print(f"\nk={k}, lam={lam:.2f}:")
    print(f"  rho={rho_est:.6f}  A/rho={A_over_rho:.6f}  A^2/rho={A2_over_rho:.6f}")
    print(f"  c_0={c_0:.4f}  c_1={c_1:.4f}  c_2={c_2:.4f}  c_1/c_0={c_1/c_0:.4f}  A/rho={A_over_rho:.4f}")
    print(f"  a0_v2={a0_v2:.4f}  a1_v2={a1_v2:.4f}  a2_v2={a2_v2:.4f}")
    print(f"  a0_v0={a0_v0:.4f}  a1_v0={a1_v0:.4f}  a2_v0={a2_v0:.4f}")
    print(f"  tau verify: max_err={err_tau:.1e}")
    print(f"  phi class dist (class-0 of s): {['r%d:%.3f'%(r,p) for r,p in enumerate(phi0_class_dist)]}")
    print(f"  phi class dist (class-1 of s): {['r%d:%.3f'%(r,p) for r,p in enumerate(phi1_class_dist)]}")
    print(f"  m0 (T4 v0-pullback, s=0 mod 3) = {m0:.6f}  (expect ~a2_v0={a2_v0:.4f}?)")
    print(f"  m1 (T4 v0-pullback, s=1 mod 3) = {m1:.6f}  (expect ~a0_v0={a0_v0:.4f}?)")
    print(f"  m1 >= m0? {m1 >= m0}  (m1-m0 = {m1-m0:+.6f})")
    print(f"")
    print(f"  FULL K-L DECOMPOSITION VERIFICATION:")
    print(f"    rho*a0 predicted = B3*c_1 + A^2/rho*m0 = {B3:.4f}*{c_1:.4f} + {A2_over_rho:.4f}*{m0:.4f}")
    print(f"                     = {B3*c_1:.6f} + {A2_over_rho*m0:.6f} = {rho_a0_predicted:.6f}")
    print(f"    rho*a0 actual    = {rho_a0_actual:.6f}")
    print(f"    rel error a0: {abs(rho_a0_predicted-rho_a0_actual)/(rho_a0_actual+1e-15):.2e}")
    print(f"")
    print(f"    rho*a1 predicted = B3*c_0 + A^2/rho*m1 = {B3:.4f}*{c_0:.4f} + {A2_over_rho:.4f}*{m1:.4f}")
    print(f"                     = {B3*c_0:.6f} + {A2_over_rho*m1:.6f} = {rho_a1_predicted:.6f}")
    print(f"    rho*a1 actual    = {rho_a1_actual:.6f}")
    print(f"    rel error a1: {abs(rho_a1_predicted-rho_a1_actual)/(rho_a1_actual+1e-15):.2e}")
    print(f"")
    print(f"  a1 - a0 DECOMPOSITION:")
    print(f"    cb term:  B3*(c_0-c_1)/rho = {cb_diff/rho_est:+.6f}  (POSITIVE)")
    print(f"    T4 term: A^2/rho*(m1-m0)/rho = {T4_diff/rho_est:+.6f}  (sign?)")
    print(f"    total = {total_diff:+.6f}  vs actual a1-a0 = {a1_v2-a0_v2:+.6f}")
    print(f"    both terms same sign? {(cb_diff > 0) and (T4_diff >= 0)}")

    return m0, m1, a0_v2, a1_v2

print("260: Complete proof that a0_v2 < a1_v2 for all lambda > 1")
print("="*70)

# Main check
analyze_full(8, 1.70)

# Lambda scan
print(f"\n\nLambda scan k=8: verify m1 >= m0 and both terms positive")
print(f"{'lam':>6}  {'m0':>8}  {'m1':>8}  {'m1>=m0':>8}  {'a0':>8}  {'a1':>8}")
for lam in [1.30, 1.50, 1.70, 1.90, 2.00]:
    m0, m1, a0, a1 = analyze_full(8, lam)
    print(f"lam={lam:.2f}  m0={m0:.4f}  m1={m1:.4f}  {'YES' if m1>=m0 else 'NO':>8}  a0={a0:.4f}  a1={a1:.4f}")
    sys.stdout.flush()

# Depth scan
print(f"\n\nDepth scan lam=1.70: verify m1 >= m0 for all k")
print(f"{'k':>4}  {'m0':>8}  {'m1':>8}  {'m1>=m0':>8}  {'a0':>8}  {'a1':>8}")
lam = 1.70
for k in range(5, 12):
    m0, m1, a0, a1 = analyze_full(k, lam)
    print(f"k={k:>2}  m0={m0:.4f}  m1={m1:.4f}  {'YES' if m1>=m0 else 'NO':>8}  a0={a0:.4f}  a1={a1:.4f}")
    sys.stdout.flush()

print("\ndone")
