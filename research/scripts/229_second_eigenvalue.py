"""
229_second_eigenvalue.py
========================
Meting f2 = (lambda2/lambda1)^2 via convergentiequotient van machtiteraten.

Methode: start power iteration vanuit v0 = ones. Na warmup:
  ratio_n = ||v_{n+1} - v_n|| / ||v_n - v_{n-1}||  ->  lambda2/lambda1

f2 = lim ratio_n^2 = spectrale convergentiefactor voor variantie.

Alternatieve meting: trac de voortgang van CV(n) = std(v_n)/mean(v_n) naar CV_inf.
  |CV(n) - CV_inf| / |CV(n-1) - CV_inf|  ->  lambda2/lambda1

Endpoint Decay: f2 < 1 <=> lambda2 < lambda1 (spectraal gat > 0).
"""
import sys
import numpy as np
from math import log2

ALPHA = log2(3.0)
LAM   = 1.70
A  = LAM ** -2.0
B1 = LAM ** (ALPHA - 2.0)
B3 = LAM ** (ALPHA - 1.0)

# Schatting CV_inf uit Script 228 (methode A, k->inf)
CV_INF_EST = 0.907   # conservatieve schatting

print(f"229: f2-meting via convergentiequotient machtiteraten  (lam={LAM})")
print(f"     A={A:.6f}  B1={B1:.6f}  B3={B3:.6f}")
print(f"     CV_inf schatting: {CV_INF_EST}")
print("=" * 65)
print()
sys.stdout.flush()

results = []

for k in range(5, 16):
    N  = 3 ** (k - 1)
    n_total = 600  # totaal iteraties
    print(f"  k={k:2d}  N={N:>10,d}  ({n_total} iter) ...", end='', flush=True)

    i  = np.arange(N, dtype=np.int64)
    T4 = (4 * i + 2) % N
    s, r = np.divmod(i, 3)
    Nl = N // 3
    m0, m2 = (r == 0), (r == 2)
    R1 = (4 * s) % Nl
    R3 = (2 * s + 1) % Nl

    def apply_F(v):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        w  = A * v[T4]
        w[m2] += B3 * cb[R3[m2]]
        w[m0] += B1 * cb[R1[m0]]
        return w

    v = np.ones(N, dtype=np.float64)
    prev_diff = None
    prev_v = None
    cv_series = []
    ratio_series = []

    for it in range(n_total):
        w = apply_F(v)
        vn = w / w.max()
        vn_norm = vn / vn.mean()  # voor CV berekening
        cv_n = float(np.std(vn_norm))
        cv_series.append(cv_n)

        if prev_v is not None:
            diff = float(np.linalg.norm(vn - prev_v))
            if prev_diff is not None and prev_diff > 1e-15 and diff > 1e-15:
                ratio = diff / prev_diff
                ratio_series.append(ratio)
            prev_diff = diff
        prev_v = vn.copy()
        v = vn

    # rho uit de laatste iteraat
    v_norm = v / v.mean()
    w_check = apply_F(v_norm)
    rho = float(w_check.mean())
    cv_final = float(np.std(v_norm))

    # Convergentiequotient: gebruik de MEEST RECENTE stabiele waarden
    ratio_arr = np.array(ratio_series)
    # Verwijder uitschieters en onbetrouwbare vroege waarden
    valid = (ratio_arr > 0.001) & (ratio_arr < 10.0)
    if valid.sum() > 100:
        # Gebruik het laatste kwart (meest uitgeconvergeerd)
        last_quarter_start = len(ratio_arr) - len(ratio_arr)//4
        recent_ratios = ratio_arr[last_quarter_start:]
        recent_valid = (recent_ratios > 0.001) & (recent_ratios < 10.0)
        if recent_valid.sum() > 10:
            ratio_mean = float(np.median(recent_ratios[recent_valid]))
        else:
            ratio_mean = float(np.median(ratio_arr[valid]))
    else:
        ratio_mean = float('nan')

    f2_est = ratio_mean ** 2

    # CV-convergentiemethode: gebruik CV_INF_EST
    cv_arr = np.array(cv_series)
    cv_gap = np.abs(cv_arr - CV_INF_EST)
    # Ratio van opeenvolgende CV-gaten: |cv_n - cv_inf| / |cv_{n-1} - cv_inf|
    cv_ratio_arr = cv_gap[1:] / (cv_gap[:-1] + 1e-30)
    # Gebruik laatste kwart
    lq = len(cv_ratio_arr) - len(cv_ratio_arr)//4
    cv_ratio_recent = cv_ratio_arr[lq:]
    cv_valid = (cv_ratio_recent > 0.01) & (cv_ratio_recent < 1.5)
    if cv_valid.sum() > 10:
        cv_ratio_median = float(np.median(cv_ratio_recent[cv_valid]))
    else:
        cv_ratio_median = float('nan')

    print(f"\n    rho={rho:.6f}  cv={cv_final:.5f}  "
          f"ratio_median={ratio_mean:.5f}  f2={f2_est:.5f}  "
          f"cv_ratio={cv_ratio_median:.5f}")

    # Extra info: std van de ratio (stabiliteit)
    if valid.sum() > 100:
        last_q = ratio_arr[last_quarter_start:]
        lq_valid = (last_q > 0.01) & (last_q < 2.0)
        print(f"    Ratio (laatste kwart): mean={float(np.mean(last_q[lq_valid])):.5f}  "
              f"median={float(np.median(last_q[lq_valid])):.5f}  "
              f"std={float(np.std(last_q[lq_valid])):.5f}  "
              f"min={float(np.min(last_q[lq_valid])):.5f}  "
              f"max={float(np.max(last_q[lq_valid])):.5f}")

    sys.stdout.flush()
    results.append({'k': k, 'N': N, 'rho': rho, 'cv': cv_final,
                    'ratio': ratio_mean, 'f2': f2_est,
                    'cv_ratio': cv_ratio_median})

print()
print()
print("=== SAMENVATTING ===")
print()
print(f"{'k':>4} {'N':>10} {'rho':>8} {'cv':>8} {'ratio':>9} {'f2':>9} {'cv_rat':>9}")
print("-" * 65)
for r in results:
    print(f"  {r['k']:>2}  {r['N']:>10,d}  {r['rho']:>8.5f}  {r['cv']:>8.5f}  "
          f"{r['ratio']:>9.5f}  {r['f2']:>9.5f}  {r['cv_ratio']:>9.5f}")

if results:
    f2_vals = [r['f2'] for r in results if not np.isnan(r['f2'])]
    cv_rat_vals = [r['cv_ratio'] for r in results if not np.isnan(r['cv_ratio'])]
    print()
    if f2_vals:
        print(f"f2 waarden:      {[f'{x:.5f}' for x in f2_vals]}")
        all_lt1 = all(x < 1.0 for x in f2_vals)
        print(f"f2 < 1 voor alle k: {'JA' if all_lt1 else 'DEELS / NEE'}")
    if cv_rat_vals:
        print(f"cv_ratio waarden: {[f'{x:.5f}' for x in cv_rat_vals]}")
        print(f"Gemiddelde cv_ratio (grote k): {np.mean(cv_rat_vals[-3:]):.5f}")

print()
print("done")
