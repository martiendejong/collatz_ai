"""
254_selfconsistent_class.py
============================
Verify the self-consistency of the mod-3 class mean structure.

FROM SCRIPT 253: Mean v2 by class:
  a0 = 0.140  (s==0 mod 3) -- SMALLEST
  a1 = 0.328  (s==1 mod 3)
  a2 = 0.363  (s==2 mod 3) -- LARGEST

HYPOTHESIS (self-consistency):
  rho * a0 ~= B3 * c1  (s==0: cb at 2*0+1=1 ==1 mod 3)
  rho * a1 ~= B3 * c0  (s==1: cb at 2*1+1=3 ==0 mod 3)
  rho * a2 ~= B3 * c2  (s==2: cb at 2*2+1=5 ==2 mod 3)

where c_r = Mean(cb[j==r mod 3]) = E[min(v2[j], v2[j+Nl3], v2[j+2Nl3]) | j==r mod 3]

KEY QUESTION: why is a0 < a1?
  a0 < a1  iff  c1 < c0  (from self-consistency)
  c_r = E[min3 from class r] < a_r always
  min-ratio k_r = c_r / a_r: higher CV => lower k_r
  CV_0 < CV_1  =>  k0 > k1  =>  c0/a0 > c1/a1

So: a0/a1 = c1/c0 -- wait no: rho*a0 = B3*c1 and rho*a1 = B3*c0
=> a0/a1 = c1/c0

And a0 < a1 iff c1 < c0. From k_r = c_r/a_r:
  c0 = k0 * a0, c1 = k1 * a1
  a0/a1 = c1/c0 = (k1*a1)/(k0*a0)
  => (a0/a1)^2 = k1/k0

So a0 < a1 iff k0 > k1 iff CV0 < CV1.

TESTS:
1. Compute c_r = Mean(cb | class r), verify rho*a_r ~= B3*c_sigma(r) where sigma = (0->1, 1->0, 2->2)
2. Compute min-ratios k_r = c_r / a_r
3. Compute CVs of v2 by class
4. Show CV0 < CV1 => k0 > k1 => self-consistent with a0 < a1
5. Compute the BETWEEN-CLASS contribution to Corr(v2[2s+1], v2[s]) and compare to actual

Also: decompose the anti-correlation into between-class and within-class parts.
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
    # Estimate rho from the K-L equation (use v at convergence)
    cb_final = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
    w_final = A * v[T4]
    w_final[m2] += B3 * cb_final[R3[m2]]
    w_final[m0] += B1 * cb_final[R1[m0]]
    # rho * v[i] = w_final[i], so rho = w_final[i]/v[i] for max element
    rho_est = float(w_final.max())
    return v, Nl, cb_final, A, B1, B3, rho_est

def analyze_selfconsistent(k, lam):
    v, Nl, cb, A, B1, B3, rho = run_kl(k, lam)
    v2 = v[2::3]
    s = np.arange(Nl, dtype=np.int64)
    Nl3 = Nl // 3

    print(f"\nk={k}, lam={lam:.2f}, Nl={Nl}, rho_est={rho:.6f}")
    print(f"  B3 = {B3:.6f}, B3/rho = {B3/rho:.6f}")

    # === CLASS MEANS AND CVs ===
    a = np.zeros(3)
    std = np.zeros(3)
    cv = np.zeros(3)
    for r in range(3):
        mask = (s % 3 == r)
        a[r] = float(np.mean(v2[mask]))
        std[r] = float(np.std(v2[mask]))
        cv[r] = std[r] / a[r]

    print(f"\n  Class means (v2):")
    for r in range(3):
        print(f"    a[{r}] = {a[r]:.6f}  std = {std[r]:.6f}  CV = {cv[r]:.4f}")
    print(f"  CV order: CV0 < CV1? {cv[0] < cv[1]}  CV0={cv[0]:.4f} CV1={cv[1]:.4f}")

    # === CB MEANS BY CLASS ===
    js = np.arange(Nl, dtype=np.int64)
    c = np.zeros(3)
    for r in range(3):
        mask = (js % 3 == r)
        c[r] = float(np.mean(cb[mask]))

    print(f"\n  cb class means:")
    for r in range(3):
        k_r = c[r] / a[r]
        print(f"    c[{r}] = {c[r]:.6f}  min-ratio k[{r}] = c[{r}]/a[{r}] = {k_r:.4f}")
    print(f"  k order: k0 > k1? {c[0]/a[0] > c[1]/a[1]}  k0={c[0]/a[0]:.4f} k1={c[1]/a[1]:.4f}")

    # === SELF-CONSISTENCY CHECK ===
    print(f"\n  Self-consistency check: rho * a[r] ~= B3 * c[sigma(r)]")
    # sigma(0)=1, sigma(1)=0, sigma(2)=2
    sigma = {0: 1, 1: 0, 2: 2}
    for r in range(3):
        lhs = rho * a[r]
        rhs = B3 * c[sigma[r]]
        print(f"    r={r}: rho*a[{r}] = {lhs:.6f}  B3*c[{sigma[r]}] = {rhs:.6f}  "
              f"ratio = {lhs/rhs:.4f}  (T4 contributes {(lhs-rhs)/lhs:.4f})")

    # === (a0/a1)^2 vs k1/k0 ===
    print(f"\n  Key ratio check: (a0/a1)^2 = {(a[0]/a[1])**2:.4f}  "
          f"k1/k0 = {(c[1]/a[1])/(c[0]/a[0]):.4f}")
    print(f"  (Should match if T4 term negligible)")

    # === BETWEEN-CLASS CONTRIBUTION TO ANTI-CORRELATION ===
    # Replace each v2[s] by its class mean a[s mod 3].
    # Then compute Corr(class_mean[2s+1 mod 3], class_mean[s mod 3]).
    max_s = (Nl - 9) // 9
    s0 = np.arange(max_s, dtype=np.int64)

    # Between-class proxy: use class means
    class_s = a[s0 % 3]
    class_2sp1 = a[(2*s0 + 1) % 3]
    corr_between = float(np.corrcoef(class_s, class_2sp1)[0,1])

    # Actual v2 anti-correlation on same set
    v2_s0 = v2[s0]
    v2_2sp1 = v2[(2*s0 + 1) % Nl]
    corr_actual = float(np.corrcoef(v2_s0, v2_2sp1)[0,1])

    print(f"\n  Anti-correlation decomposition:")
    print(f"    Actual Corr(v2[2s+1], v2[s])           = {corr_actual:+.4f}")
    print(f"    Between-class Corr (using class means)  = {corr_between:+.4f}")
    print(f"    Fraction explained by between-class     = {corr_between/corr_actual:.4f}")

    # Note: the between-class corr using discrete {a0, a1, a2} is:
    # Corr_btwn = (1/3*[a0*a1 + a1*a0 + a2*a2] - mean_a^2) / Var_btwn
    # where mean_a = (a0+a1+a2)/3 and Var_btwn = (a0^2+a1^2+a2^2)/3 - mean_a^2
    mean_a = (a[0]+a[1]+a[2])/3
    var_btwn = (a[0]**2 + a[1]**2 + a[2]**2)/3 - mean_a**2
    cov_btwn = (a[0]*a[1] + a[1]*a[0] + a[2]*a[2])/3 - mean_a**2
    corr_btwn_theory = cov_btwn / var_btwn
    print(f"    Between-class Corr (analytic formula)   = {corr_btwn_theory:+.4f}")
    print(f"    (Formula: [2*a0*a1 + a2^2]/3 - mean_a^2 / Var(class means))")

    # Contribution from each pair:
    pairs = [(0,1,'0->1'), (1,0,'1->0'), (2,2,'2->2')]
    print(f"\n  Pair-by-pair: Corr(v2[2s+1], v2[s]) split by s mod 3:")
    for r_s, r_2sp1, label in pairs:
        mask = (s0 % 3 == r_s)
        if mask.sum() > 2:
            c_pair = float(np.corrcoef(v2_s0[mask], v2_2sp1[mask])[0,1])
            between_effect = (a[r_2sp1] - a[r_s]) / (np.sqrt(np.var(v2_s0[mask])) * np.sqrt(np.var(v2_2sp1[mask])))
            print(f"    {label}: Corr = {c_pair:+.4f}  "
                  f"(mean shift: {a[r_s]:.3f} -> {a[r_2sp1]:.3f})  n={mask.sum()}")

    return corr_actual, corr_between, corr_btwn_theory

print("254: Self-consistency of mod-3 class structure")
print("="*70)

analyze_selfconsistent(8, 1.70)

print("\n\n=== Lambda scan k=8: self-consistency and anti-corr decomposition ===")
print(f"{'lam':>6}  {'a0':>7}  {'a1':>7}  {'a0/a1':>7}  "
      f"{'k0':>7}  {'k1':>7}  {'k1/k0':>7}  {'(a0/a1)^2':>10}  {'corr_act':>9}  {'corr_btwn':>10}")
for lam in [1.30, 1.50, 1.70, 1.90, 2.00]:
    v, Nl, cb, A, B1, B3, rho = run_kl(8, lam)
    v2 = v[2::3]
    js = s = np.arange(Nl, dtype=np.int64)
    a = np.array([float(np.mean(v2[s%3==r])) for r in range(3)])
    c = np.array([float(np.mean(cb[s%3==r])) for r in range(3)])
    k_r = c / a  # min-ratio by class
    max_s = (Nl-9)//9
    s0 = np.arange(max_s, dtype=np.int64)
    c_actual = float(np.corrcoef(v2[s0], v2[(2*s0+1)%Nl])[0,1])
    mean_a = a.mean(); var_btwn = np.mean(a**2) - mean_a**2
    cov_btwn = (a[0]*a[1] + a[1]*a[0] + a[2]*a[2])/3 - mean_a**2
    c_btwn = cov_btwn / var_btwn
    print(f"lam={lam:.2f}  {a[0]:>7.4f}  {a[1]:>7.4f}  {a[0]/a[1]:>7.4f}  "
          f"{k_r[0]:>7.4f}  {k_r[1]:>7.4f}  {k_r[1]/k_r[0]:>7.4f}  "
          f"{(a[0]/a[1])**2:>10.4f}  {c_actual:>9.4f}  {c_btwn:>10.4f}")
    sys.stdout.flush()

print("\ndone")
