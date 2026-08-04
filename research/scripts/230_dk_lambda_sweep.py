"""
230_dk_lambda_sweep.py
======================
Meet d_k(lambda) voor lambda in {1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.5, 3.0}
en k=13..17.

Doel: fit d_inf(lambda) als functie van A, B1, B3 en zoek een analytische
formule. Obs 401 stelde separabiliteit d_k(lam) ~ f(lam) + g(k) voor.

d_k = var_end(k+1)/var_end(k).
var_end(k) = Var_s[log2(T[:,s]) - log2(mean(T[:,s]))] met T[:,s]=[v[s],v[s+Nl],v[s+2Nl]].
"""
import sys
import numpy as np
from math import log2

ALPHA = log2(3.0)  # log2(3) ≈ 1.585

LAMBDAS = [1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.5, 3.0]
K_RANGE = list(range(13, 18))  # k=13..17
N_ITER  = 300

print("230: d_k(lambda) sweep")
print(f"     lambda = {LAMBDAS}")
print(f"     k      = {K_RANGE}")
print("=" * 70)
print()
sys.stdout.flush()


def var_end_from_v(v, k):
    Nl = 3 ** (k - 2)
    T = np.stack([v[:Nl], v[Nl:2*Nl], v[2*Nl:]])
    lmean = T.mean(axis=0)
    X = np.log2(T) - np.log2(lmean)[None, :]
    return float(np.var(X))


def compute_eigvec(k, lam, n_iter=N_ITER):
    A  = lam ** -2.0
    B1 = lam ** (ALPHA - 2.0)
    B3 = lam ** (ALPHA - 1.0)
    N  = 3 ** (k - 1)
    Nl = N // 3
    i  = np.arange(N, dtype=np.int64)
    T4 = (4 * i + 2) % N
    s, r = np.divmod(i, 3)
    m0, m2 = (r == 0), (r == 2)
    R1 = (4 * s) % Nl
    R3 = (2 * s + 1) % Nl
    v  = np.ones(N, dtype=np.float64)
    for _ in range(n_iter):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        w  = A * v[T4]
        w[m2] += B3 * cb[R3[m2]]
        w[m0] += B1 * cb[R1[m0]]
        v = w / w.max()
    return v


# Collect var_end(k, lam) for all k and lam
# Then d_k(lam) = var_end(k+1,lam)/var_end(k,lam)

results = {}  # (lam, k) -> var_end

for lam in LAMBDAS:
    A  = lam ** -2.0
    B1 = lam ** (ALPHA - 2.0)
    B3 = lam ** (ALPHA - 1.0)
    print(f"lambda={lam:.2f}  A={A:.5f}  B1={B1:.5f}  B3={B3:.5f}")
    sys.stdout.flush()

    for k in K_RANGE:
        v  = compute_eigvec(k, lam)
        ve = var_end_from_v(v, k)
        results[(lam, k)] = ve
        print(f"  k={k}: var_end={ve:.8f}")
        sys.stdout.flush()

    print()

# Bereken d_k(lam) = var_end(k+1)/var_end(k) voor k=13..16
print()
print("=" * 70)
print("d_k(lambda) TABEL")
print()

# Header
header = f"{'lam':>6}  {'A':>7}  {'B1':>7}  {'B3':>7}"
for k in K_RANGE[:-1]:
    header += f"  {'d_'+str(k):>8}"
header += f"  {'d_avg':>8}"
print(header)
print("-" * len(header))

dk_table = {}  # (lam, k) -> d_k

for lam in LAMBDAS:
    A  = lam ** -2.0
    B1 = lam ** (ALPHA - 2.0)
    B3 = lam ** (ALPHA - 1.0)
    row = f"{lam:6.2f}  {A:7.5f}  {B1:7.5f}  {B3:7.5f}"
    dk_vals = []
    for k in K_RANGE[:-1]:
        ve_k  = results[(lam, k)]
        ve_k1 = results[(lam, k+1)]
        dk    = ve_k1 / ve_k
        dk_table[(lam, k)] = dk
        row  += f"  {dk:8.6f}"
        dk_vals.append(dk)
    avg = float(np.mean(dk_vals))
    row += f"  {avg:8.6f}"
    print(row)

print()

# Extrapolatie: d_inf(lam) via incrementenreeks d_k - d_{k-1}
print("=" * 70)
print("d_inf(lambda) EXTRAPOLATIE (lineaire fit op d_k vs k)")
print()
print(f"{'lam':>6}  {'d_inf_lin':>10}  {'d_15':>8}  {'d_16':>8}  {'alle<1':>7}")
print("-" * 50)

for lam in LAMBDAS:
    ks = K_RANGE[:-1]  # k=13..16
    dk_vals = [dk_table[(lam, k)] for k in ks]
    # Lineaire fit d_k = a + b*k
    coeffs = np.polyfit(ks, dk_vals, 1)
    slope, intercept = float(coeffs[0]), float(coeffs[1])
    d_inf_lin = intercept + slope * 100  # extrapolatie naar k=100 (ruwe schatting)
    d15 = dk_table[(lam, 15)]
    d16 = dk_table[(lam, 16)]
    alle_lt1 = 'JA' if all(d < 1.0 for d in dk_vals) else 'NEE'
    print(f"{lam:6.2f}  {d_inf_lin:10.6f}  {d15:8.6f}  {d16:8.6f}  {alle_lt1:>7}")

print()

# Fit d_inf(lam) als functie van A, B1, B3
# Probeer: d_inf = 1 - c * (B1 + B3) / (A + B1 + B3)  [heuristiek]
print("=" * 70)
print("d_gemiddeld(lam) vs lambda-eigenschappen (heuristiek fit)")
print()
print(f"{'lam':>6}  {'d_avg':>8}  {'B_frac':>8}  {'A_lam':>8}  {'1-B_frac':>10}")
print("-" * 55)
for lam in LAMBDAS:
    A  = lam ** -2.0
    B1 = lam ** (ALPHA - 2.0)
    B3 = lam ** (ALPHA - 1.0)
    ks = K_RANGE[:-1]
    dk_vals = [dk_table[(lam, k)] for k in ks]
    d_avg = float(np.mean(dk_vals))
    B_tot = B1 + B3
    A_tot = A + B1 + B3
    B_frac = B_tot / A_tot
    print(f"{lam:6.2f}  {d_avg:8.6f}  {B_frac:8.5f}  {A:8.5f}  {1-B_frac:10.5f}")

print()
print("done")
