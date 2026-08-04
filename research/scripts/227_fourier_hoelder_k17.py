"""
227_fourier_hoelder_k17.py
===========================
Fourier-Hoelder-exponent van de K-L Perron-eigenvector voor k=14..17.

Verlenging van Obs 424 (k=10,12,14,15,16).
Doel: aanscherpen van alpha_inf >= 0.646 (Conjecture G numeriek).

Eerdere resultaten:
  k=10: |alpha| = 0.706
  k=12: |alpha| = 0.687
  k=14: |alpha| = 0.675
  k=15: |alpha| = 0.670
  k=16: |alpha| = 0.666

VECTORIZED v3-berekening: geen pure Python loop, elk niveau via numpy.
"""
import sys
import numpy as np
from math import log2, log

ALPHA_MATH = log2(3.0)
LAM   = 1.70
A  = LAM ** -2.0
B1 = LAM ** (ALPHA_MATH - 2.0)
B3 = LAM ** (ALPHA_MATH - 1.0)


def compute_v3_vectorized(N):
    """Bereken 3-adische valuatie voor n=0..N-1 via numpy (O(k) passes)."""
    v3 = np.zeros(N, dtype=np.int32)
    m = np.arange(N, dtype=np.int64)
    mask = (m > 0)   # skip n=0
    while True:
        div_mask = mask & (m % 3 == 0)
        if not div_mask.any():
            break
        m[div_mask] //= 3
        v3[div_mask] += 1
        mask = div_mask  # alleen doorwerken op al-deelbare
    return v3   # v3[0] = 999 nee, = 0 maar n=0 skippen we via masker


def compute_alpha_k(k, n_iter=None):
    """Bereken Hoelder-exponent alpha voor grootte k."""
    N  = 3 ** (k - 1)
    if n_iter is None:
        n_iter = 300 if k <= 15 else 200  # minder iteraties voor grote k

    i  = np.arange(N, dtype=np.int64)
    T4 = (4 * i + 2) % N
    s, r = np.divmod(i, 3)
    Nl = N // 3
    m0, m2 = (r == 0), (r == 2)
    R1 = (4 * s) % Nl
    R3 = (2 * s + 1) % Nl

    print(f"  k={k}: N={N:,d}, iter={n_iter} ...")
    sys.stdout.flush()

    # Perron-eigenvector
    v = np.ones(N, dtype=np.float64)
    for it in range(n_iter):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        w  = A * v[T4]
        w[m2] += B3 * cb[R3[m2]]
        w[m0] += B1 * cb[R1[m0]]
        v = w / w.max()
        if it % 100 == 99:
            print(f"    iter {it+1}/{n_iter} ...")
            sys.stdout.flush()
    v /= v.mean()
    print(f"  k={k}: eigenvector gereed, CV={float(np.std(v)):.5f}")
    sys.stdout.flush()

    # FFT
    print(f"  k={k}: FFT berekenen ...")
    sys.stdout.flush()
    vhat = np.fft.rfft(v)   # rfft is sneller dan fft voor reele input
    abs_vhat = np.abs(vhat)
    # rfft geeft indices 0..N//2; we gebruiken 1..N//2 (skip dc-component)
    n_rfft = len(abs_vhat)
    print(f"  k={k}: FFT klaar ({n_rfft} coefficienten)")
    sys.stdout.flush()

    # Bereken v3 voor n=1..n_rfft-1 vectorized
    print(f"  k={k}: v3-berekening ...")
    sys.stdout.flush()
    ns = np.arange(n_rfft, dtype=np.int64)
    # v3 berekenen:
    v3_arr = np.zeros(n_rfft, dtype=np.int32)
    temp = ns.copy()
    active = (temp > 0)
    while active.any():
        new_active = active & (temp % 3 == 0)
        if not new_active.any():
            break
        temp[new_active] //= 3
        v3_arr[new_active] += 1
        active = new_active
    print(f"  k={k}: v3-berekening klaar, max v3={v3_arr[1:].max()}")
    sys.stdout.flush()

    # Groepeer per v3-niveau (skip n=0 want n=0 is dc)
    levels = {}
    for lev in range(k):
        mask_lev = (ns > 0) & (v3_arr == lev)
        if mask_lev.any():
            log_vhat_lev = np.log(abs_vhat[mask_lev] + 1e-300)
            levels[lev] = float(np.mean(log_vhat_lev))

    pts = []
    for lev in sorted(levels.keys()):
        log_norm = -lev * log(3.0)  # log(3^{-lev}) = -lev * log3
        pts.append((log_norm, levels[lev]))

    log_norms = np.array([p[0] for p in pts])
    log_vhats = np.array([p[1] for p in pts])

    ok = np.isfinite(log_vhats) & (log_norms > -15)
    if ok.sum() >= 3:
        slope, intercept = np.polyfit(log_norms[ok], log_vhats[ok], 1)
        alpha = slope
    else:
        alpha = float('nan')

    print(f"  k={k}: alpha = {alpha:.4f} (via rfft, {len(pts)} niveaus)")
    print(f"  Niveaus (log_norm, log_vhat): {[(round(p[0],2), round(p[1],3)) for p in pts[:8]]}")
    sys.stdout.flush()
    return alpha


print(f"227: Fourier-Hoelder-exponent k=14..17  (lam={LAM})")
print("=" * 60)
print()
print("Eerdere resultaten (Obs 424):")
print("  k=10: 0.706, k=12: 0.687, k=14: 0.675, k=15: 0.670, k=16: 0.666")
print("  Decrements: 0.019(10->12), 0.012(12->14), 0.005(14->15), 0.004(15->16)")
print("  Ratio decrements: ~0.8; extrapolatie alpha_inf >= 0.646")
print()
sys.stdout.flush()

results = {}
for k in range(14, 18):
    alpha = compute_alpha_k(k)
    results[k] = alpha
    print()
    sys.stdout.flush()

print("=== SAMENVATTING ===")
ks = sorted(results.keys())
alphas = [results[k] for k in ks]
print(f"{'k':>4}  {'alpha':>8}")
for k, a in zip(ks, alphas):
    print(f"  {k:2d}  {a:.4f}")

print()
decrements = np.diff(alphas)
print(f"Decrements (nieuw): {[f'{d:.5f}' for d in decrements]}")
if len(decrements) >= 2:
    dec_ratios = decrements[1:] / (np.abs(decrements[:-1]) + 1e-12)
    print(f"Ratio decrements:    {[f'{r:.4f}' for r in dec_ratios]}")

# Gecombineerd met eerdere data
all_ks    = [10, 12, 14, 15, 16] + [k for k in ks if k > 16]
all_alpha = [0.706, 0.687, 0.675, 0.670, 0.666] + [results[k] for k in ks if k > 16]
print()
print("Gecombineerde serie:")
prev_a = None
for k2, a2 in zip(all_ks, all_alpha):
    dcv = f"  dec={a2-prev_a:+.4f}" if prev_a is not None else ""
    print(f"  k={k2:2d}  alpha={a2:.4f}{dcv}")
    prev_a = a2

# Extrapolatie verfijnd
if len(all_alpha) >= 2:
    all_decs = np.diff(all_alpha)
    print()
    print(f"Alle decrements: {[f'{d:.4f}' for d in all_decs]}")
    all_dec_ratios = all_decs[1:] / (np.abs(all_decs[:-1]) + 1e-12)
    print(f"Ratio decrements: {[f'{r:.4f}' for r in all_dec_ratios]}")
    avg_ratio = float(np.mean(all_dec_ratios[-3:]))
    last_dec = abs(all_decs[-1])
    if avg_ratio < 1.0:
        tail = last_dec / (1.0 - avg_ratio)
        alpha_inf_lb = all_alpha[-1] - tail
        print(f"Extrapolatie (ratio={avg_ratio:.4f}): alpha_inf >= {all_alpha[-1]:.4f} - {tail:.4f} = {alpha_inf_lb:.4f}")
    print()
    print(f"SCHATTING: alpha_inf ~ {all_alpha[-1]:.4f} (bovengrens)")
    print(f"           alpha_inf >= {all_alpha[-1] - last_dec:.4f} (ondergrens, 1 decrement)")

sys.stdout.flush()
print()
print("done")
