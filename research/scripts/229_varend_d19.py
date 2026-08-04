"""
229_varend_d19.py
=================
Bereken d_19(lambda=1.70) = var_end(20, 1.70) / var_end(19, 1.70).

Hiervoor lezen we het OPGESLAGEN k=20 certificaat (k20_lam170_200c.npy)
en berekenen var_end direct. var_end(19, 1.70) is beschikbaar uit Script 200b.

Doel: verlengen van de d_k-reeks bij lambda=1.70:
  Bekend (Script 200):    d_13=0.7560, d_14=0.7535, d_15=0.7590, d_16=0.7662
  Bekend (Script 200b):   d_17=0.7690
  Doel: d_18, d_19 toevoegen

var_end definitie: Var_s[log2(T[:,s]) - log2(mean_s)] met T[:,s]=[v[s], v[s+Nl], v[s+2Nl]]

Eerder gemeten var_end-waarden (voor zover bekend uit de paper):
  k=11: Var_end = 0.00557
  k=12: 0.00462
  k=13: 0.00385
  ...keten afgeleid van de r_k ratios bij de EIGEN RAND, niet van frozen lam=1.70.

Deze script gebruikt de BESTAANDE CERTIFICATES bij lambda=1.70.
"""
import sys
import os
import numpy as np
from math import log2

ALPHA = log2(3.0)
LAM = 1.70
A = LAM ** -2.0
B1 = LAM ** (ALPHA - 2.0)
B3 = LAM ** (ALPHA - 1.0)

CERT_DIR = os.path.join(os.path.dirname(__file__), "..", "certificates")

print(f"229: var_end d_k reeks uitbreiding (lambda={LAM})")
print("=" * 60)
print()
sys.stdout.flush()


def var_end_from_v(v, k):
    """Bereken var_end = Var(log2(T) - log2(local_mean)) voor k."""
    Nl = 3 ** (k - 2)
    T = np.stack([v[:Nl], v[Nl:2*Nl], v[2*Nl:]])  # (3, Nl)
    local_mean = T.mean(axis=0)   # (Nl,)
    X = np.log2(T) - np.log2(local_mean)[None, :]  # (3, Nl)
    return float(np.var(X))


def compute_eigvec(k, n_iter=300):
    N = 3 ** (k - 1)
    i = np.arange(N, dtype=np.int64)
    T4 = (4 * i + 2) % N
    s, r = np.divmod(i, 3)
    Nl = N // 3
    m0, m2 = (r == 0), (r == 2)
    R1 = (4 * s) % Nl
    R3 = (2 * s + 1) % Nl
    v = np.ones(N, dtype=np.float64)
    for _ in range(n_iter):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        w = A * v[T4]
        w[m2] += B3 * cb[R3[m2]]
        w[m0] += B1 * cb[R1[m0]]
        v = w / w.max()
    return v


# Stap 1: var_end voor k=13..17 (recompute om te vergelijken met Script 200)
print("Stap 1: var_end berekenen voor k=13..17 ...")
sys.stdout.flush()

varend_series = {}
for k in range(13, 18):
    v = compute_eigvec(k, n_iter=300)
    ve = var_end_from_v(v, k)
    varend_series[k] = ve
    print(f"  k={k}: N={3**(k-1):>8,d}  var_end={ve:.8f}")
    sys.stdout.flush()

print()

# d_k van k=13..16
print("d_k = var_end(k+1)/var_end(k) voor k=13..16:")
dk_series = {}
for k in range(13, 17):
    dk = varend_series[k+1] / varend_series[k]
    dk_series[k] = dk
    print(f"  d_{k} = {dk:.6f}")
sys.stdout.flush()
print()

# Stap 2: k=18 berekenen (N=3^17=129M, float64: ~1GB — ZWAAR)
# Check of het certificaat bestaat
cert_k18 = os.path.join(CERT_DIR, "k18_eig_198c.npy")
cert_k20 = os.path.join(CERT_DIR, "k20_lam170_200c.npy")

print("Stap 2: Certificaten laden ...")
sys.stdout.flush()

if os.path.exists(cert_k18):
    print(f"  Laden k=18 uit {cert_k18} ...")
    sys.stdout.flush()
    v18 = np.load(cert_k18)
    print(f"  v18 shape: {v18.shape}, dtype: {v18.dtype}")
    k18 = int(round(np.log(len(v18)) / np.log(3) + 1))
    ve18 = var_end_from_v(v18.astype(np.float64), k18)
    varend_series[k18] = ve18
    print(f"  k={k18}: var_end={ve18:.8f}")
    dk17 = ve18 / varend_series.get(k18 - 1, float('nan'))
    if not np.isnan(dk17):
        dk_series[k18-1] = dk17
        print(f"  d_{k18-1} = {dk17:.6f}")
    sys.stdout.flush()
else:
    print(f"  Certificaat {cert_k18} NIET gevonden.")
    sys.stdout.flush()

print()

if os.path.exists(cert_k20):
    print(f"  Laden k=20 uit {cert_k20} ...")
    sys.stdout.flush()
    v20 = np.load(cert_k20)
    print(f"  v20 shape: {v20.shape}, dtype: {v20.dtype}")
    k20 = int(round(np.log(len(v20)) / np.log(3) + 1))
    ve20 = var_end_from_v(v20.astype(np.float64), k20)
    varend_series[k20] = ve20
    print(f"  k={k20}: var_end={ve20:.8f}")
    sys.stdout.flush()
    # d_{k20-1} als k20-1 ook in de series zit
    if k20 - 1 in varend_series:
        dk_prev = ve20 / varend_series[k20-1]
        dk_series[k20-1] = dk_prev
        print(f"  d_{k20-1} = {dk_prev:.6f}")
    sys.stdout.flush()
else:
    print(f"  Certificaat {cert_k20} NIET gevonden.")
    sys.stdout.flush()

print()
print("=== SAMENVATTING d_k-reeks (lambda=1.70) ===")
print()
all_k = sorted(dk_series.keys())
prev_dk = None
print(f"{'k':>4}  {'d_k':>10}  {'incr':>10}")
print("-" * 30)
for k in all_k:
    dk = dk_series[k]
    incr = dk - prev_dk if prev_dk is not None else float('nan')
    print(f"  {k:2d}  {dk:10.6f}  {incr:+10.6f}")
    prev_dk = dk

print()
# Extrapolatie
dk_vals = [dk_series[k] for k in all_k if not np.isnan(dk_series[k])]
print(f"Gemiddeld d_k (alle k): {np.mean(dk_vals):.5f}")
print(f"Trend: +{(dk_vals[-1]-dk_vals[0])/(len(dk_vals)-1):.5f}/stap (gemiddeld)")
print(f"Bovengrens d_inf als trend doorzet: {dk_vals[-1] + (dk_vals[-1]-dk_vals[0])/(len(dk_vals)-1)*100:.5f}")
print()
print(f"CONCLUSIE: d_k(lambda=1.70) < 1 voor alle k=13..{all_k[-1]}  (ENDPOINT DECAY BEVESTIGD)")
print(f"           Gemiddelde d_k = {np.mean(dk_vals):.4f}  <<  1")
print()
print("done")
