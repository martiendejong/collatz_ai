"""
131_orbit_length_dist.py
=========================
Distribution of Collatz orbit lengths (in macro-steps).

From the Lyapunov analysis (Obs 268):
- Mean log-ratio per step: mu = -0.575
- Std per step: sigma = 11.1
- Orbit length T ~ Total_decrease / |mu| = b*log(2) / 0.575 ~ 1.21*b
- By CLT, T should be approximately Gaussian with mean 1.21b and std?

Key question: Is the orbit length distribution actually Gaussian?
What's the variance of T?

CLT for stopping times:
If X_t are i.i.d. with mean mu<0 and variance sig2, and T = first t where
sum_{i=1}^T X_i <= -L, then by the central limit theorem for random walks:
  T ~ L/|mu| (linear in L=b*log2)
  Var(T) ~ L * sig2 / mu^3 (?) -- need Wald's identity for random walks

Actually by Wald's identity: E[S_T] = E[T] * mu (where S_T = sum of X_i up to T)
And Var(S_T) = E[T] * sig2 + Var(T) * mu^2 (by Wald's second identity for stopped walks)
But S_T ≈ -L (stopping criterion), so Var(S_T) ≈ 0 (approximately), giving:
  Var(T) ≈ E[T] * sig2 / mu^2

std(T) ≈ sqrt(E[T] * sig2 / mu^2) = sqrt(T * sig2) / |mu|
       = sqrt(1.21*b * 11.1^2) / 0.575
       = 11.1 * sqrt(1.21*b) / 0.575
       ~ 21.3 * sqrt(b)

So the orbit length distribution should be:
  T ~ Normal(1.21*b, (21.3*sqrt(b))^2)
  std(T) ~ 21.3 * sqrt(b)

This is a key prediction to verify.
"""
import numpy as np
import math
from scipy import stats
import random as _r

def v2(x):
    if x == 0: return 999
    c = 0
    while x % 2 == 0: x //= 2; c += 1
    return c

def macro_step(n):
    K = v2(n+1); m = (n+1) >> K; x = m*(3**K)-1; l = v2(x)
    return x >> l, K, l

def orbit_length(n):
    """Number of macro-steps until n reaches 1."""
    t = 0
    while n > 1:
        n, _, _ = macro_step(n)
        t += 1
    return t

print("=" * 70)
print("PART 1: ORBIT LENGTH DISTRIBUTION FOR VARIOUS BIT LENGTHS")
print("=" * 70)
print()

_r.seed(42)
MU = 2*math.log(3) - 4*math.log(2)  # -0.5754
SIGMA = 11.12  # measured in script 124
B_LN2 = math.log(2)  # log(2)

results = {}
for b in [30, 50, 100, 200, 500]:
    N_samples = 500
    lengths = []
    for _ in range(N_samples):
        n = _r.getrandbits(b) | 1  # random b-bit odd
        lengths.append(orbit_length(n))
    T = np.array(lengths)
    results[b] = T

    T_mean = b * B_LN2 / abs(MU)  # predicted mean
    T_std = SIGMA * math.sqrt(T_mean) / abs(MU)  # predicted std
    print(f"b={b:>4}: n={N_samples} orbits")
    print(f"  Observed:  mean={T.mean():.1f}, std={T.std():.1f}")
    print(f"  Predicted: mean={T_mean:.1f}, std={T_std:.1f}")
    print(f"  (mean/b)={T.mean()/b:.4f} vs {T_mean/b:.4f}")
    print(f"  (std/sqrt(b))={T.std()/math.sqrt(b):.3f} vs {T_std/math.sqrt(b):.3f}")
    print()

print("=" * 70)
print("PART 2: NORMALITY TEST OF ORBIT LENGTH DISTRIBUTION")
print("=" * 70)
print()

for b in [50, 100, 200]:
    T = results[b]
    T_mean = b * B_LN2 / abs(MU)
    T_std = SIGMA * math.sqrt(T_mean) / abs(MU)
    # Shapiro-Wilk test for normality
    stat, pval = stats.shapiro(T)
    # Kolmogorov-Smirnov test against Normal
    T_standardized = (T - T.mean()) / T.std()
    ks_stat, ks_pval = stats.kstest(T_standardized, 'norm')
    # Skewness and kurtosis
    skew = stats.skew(T)
    kurt = stats.kurtosis(T)  # excess kurtosis (0 for normal)
    print(f"b={b}: Shapiro-Wilk p={pval:.4f}, KS p={ks_pval:.4f}")
    print(f"  Skewness={skew:.4f} (0=normal)")
    print(f"  Excess kurtosis={kurt:.4f} (0=normal)")
    print()

print("=" * 70)
print("PART 3: ORBIT LENGTH SCALING FIT")
print("=" * 70)
print()

b_list = [30, 50, 100, 200, 500]
mean_T = [results[b].mean() for b in b_list]
std_T = [results[b].std() for b in b_list]

# Fit mean: T_mean = alpha * b + beta
log_b = np.log(b_list)
log_mean = np.log(mean_T)
slope_mean, intercept_mean = np.polyfit(log_b, log_mean, 1)
print(f"Mean orbit length scaling: T_mean ~ b^{slope_mean:.4f}")
print(f"  Coefficient: {math.exp(intercept_mean):.4f}")
print(f"  Expected: T_mean = {B_LN2/abs(MU):.4f} * b = {B_LN2/abs(MU):.4f}*b")

# Fit std: T_std = c * b^gamma
log_std = np.log(std_T)
slope_std, intercept_std = np.polyfit(log_b, log_std, 1)
print(f"\nStd orbit length scaling: T_std ~ b^{slope_std:.4f}")
print(f"  Coefficient: {math.exp(intercept_std):.4f}")
T_mean_avg = np.mean(mean_T) / np.mean(b_list)
expected_std_coeff = SIGMA / abs(MU) / math.sqrt(B_LN2/abs(MU))
print(f"  Expected: T_std ~ {SIGMA/abs(MU):.3f} * sqrt(T_mean) = {SIGMA/abs(MU):.3f} * sqrt({B_LN2/abs(MU):.3f}*b)")
print(f"           = {SIGMA/abs(MU) * math.sqrt(B_LN2/abs(MU)):.3f} * sqrt(b)")

print()
print("=" * 70)
print("PART 4: ORBIT LENGTH QUANTILES vs GAUSSIAN PREDICTION")
print("=" * 70)
print()

for b in [100, 200]:
    T = results[b]
    T_mean_pred = b * B_LN2 / abs(MU)
    T_std_pred = SIGMA * math.sqrt(T_mean_pred) / abs(MU)
    print(f"b={b}: Predicted N({T_mean_pred:.0f}, {T_std_pred:.0f}^2):")
    print(f"  {'Quantile':>10} {'Predicted':>12} {'Observed':>12} {'Diff':>8}")
    for q in [0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95]:
        pred = T_mean_pred + T_std_pred * stats.norm.ppf(q)
        obs = np.quantile(T, q)
        print(f"  {q:>10.2f} {pred:>12.1f} {obs:>12.1f} {obs-pred:>8.1f}")
    print()

print("=" * 70)
print("PART 5: VARIANCE FORMULA VERIFICATION")
print("=" * 70)
print()
print("Wald's second moment identity for random walk stopped at -L:")
print("  Var(T) = E[T] * Var(X) / (E[X])^2")
print("         = (b*log2/|mu|) * sigma^2 / mu^2")
print(f"         = (b*{B_LN2:.4f}/{abs(MU):.4f}) * {SIGMA**2:.2f} / {MU**2:.4f}")
print(f"         = b * {B_LN2/abs(MU) * SIGMA**2 / MU**2:.2f}")
print(f"  Std(T) = sqrt(b) * {math.sqrt(B_LN2/abs(MU) * SIGMA**2 / MU**2):.2f}")
print()
print("Observed std vs predicted:")
for b in b_list:
    T = results[b]
    pred_var = b * B_LN2/abs(MU) * SIGMA**2 / MU**2
    pred_std = math.sqrt(pred_var)
    print(f"  b={b:>4}: pred_std={pred_std:.1f}, obs_std={T.std():.1f}, ratio={T.std()/pred_std:.3f}")
