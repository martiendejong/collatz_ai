"""
232_sigma1_autocorr.py
======================
Meet de sigma1-autocorrelatie van log v(r=0) voor k=13..16 bij lambda=1.70.

sigma1(s) = (4s+2) mod Nl is een enkelvoudige Nl-cyclus (Lem lem:sigma1).
De "mixing factor" is:
  rho1(k) = Cov_s[f0(s), f0(sigma1(s))] / Var_s[f0(s)]
           waar f0(s) = log2 v(r=0, s)

Als rho1 < 1: de sigma1-shuffle verlaagt de correlatie tussen f0 op
aangrenzende posities in de sigma1-cyclische ordening. Dit is de kwantitatieve
grondslag voor d_k < 1.

Aanvullende metingen:
  rho_cross(k) = Cov_s[f0(s), f1(s)] / sqrt(Var[f0]*Var[f1])
               = Cov[f0(s), f0(sigma1(s))] / sqrt(Var[f0]*Var[f0])  [want f1=f0 circ sigma1 + const]
               = rho1  (check: moeten gelijk zijn)

  delta_f = Var[f0(sigma1(s)) - f0(s)] / Var[f0(s)] = 2*(1 - rho1)  (mixing gain)

  Decompose var_end in r=0, r=1, r=2 bijdragen:
  ve0 = Var_s[X0(s)], ve1 = Var_s[X1(s)], ve2 = Var_s[X2(s)]
  var_end = (ve0 + ve1 + ve2)/3  [bij benadering, als mean(X_r)~0]

Bewijsrelevantie:
  Als rho1 < 1 (en constant met k), dan biedt de sigma1-single-cycle een
  ANALYTISCH BEWIJS van d_k < 1 via Cauchy-Schwarz + ergodische mixing.
"""
import sys
import numpy as np
from math import log2

ALPHA = log2(3.0)
LAM   = 1.70
A  = LAM ** -2.0
B1 = LAM ** (ALPHA - 2.0)
B3 = LAM ** (ALPHA - 1.0)
N_ITER = 300

print("232: sigma1-autocorrelatie en mixing factor (lambda=1.70)")
print(f"     A={A:.6f}  B1={B1:.6f}  B3={B3:.6f}  N_iter={N_ITER}")
print("=" * 70)
print()
sys.stdout.flush()

K_RANGE = [13, 14, 15, 16]

results = {}

for k in K_RANGE:
    N  = 3 ** (k - 1)
    Nl = N // 3
    print(f"k={k}: N={N:,}  Nl={Nl:,}")
    sys.stdout.flush()

    # Indexen
    i  = np.arange(N, dtype=np.int64)
    T4 = (4 * i + 2) % N
    s, r = np.divmod(i, 3)
    m0, m2 = (r == 0), (r == 2)
    R1 = (4 * s) % Nl
    R3 = (2 * s + 1) % Nl

    # sigma1 op {0,...,Nl-1}
    sl = np.arange(Nl, dtype=np.int64)
    sigma1 = (4 * sl + 2) % Nl   # sigma1(s)

    # Power-iteratie
    v = np.ones(N, dtype=np.float64)
    for it in range(N_ITER):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        w  = A * v[T4]
        w[m2] += B3 * cb[R3[m2]]
        w[m0] += B1 * cb[R1[m0]]
        v = w / w.max()

    # Eigenvector verdeling over de drie types
    v0 = v[:Nl]       # r=0
    v1 = v[Nl:2*Nl]   # r=1
    v2 = v[2*Nl:]     # r=2

    # log-waarden
    f0 = np.log2(v0)
    f1 = np.log2(v1)
    f2 = np.log2(v2)

    # Lokaal gemiddelde (rekenkundig van de drie types)
    m_s = (v0 + v1 + v2) / 3.0    # (Nl,) lokaal gemiddelde
    lm  = np.log2(m_s)

    # Within-triplet log-deviaties: X_r(s) = log2(v_r(s)) - log2(m(s))
    X0 = f0 - lm
    X1 = f1 - lm
    X2 = f2 - lm

    var_end = float(np.var(np.stack([X0, X1, X2])))
    ve0 = float(np.var(X0))
    ve1 = float(np.var(X1))
    ve2 = float(np.var(X2))

    # --- Sigma1-autocorrelatie ---
    # f0(sigma1(s)): haalt v(r=0) op bij positie sigma1(s)
    f0_shifted = f0[sigma1]   # f0 geëvalueerd op sigma1(s)

    # Autocorrelatie bij lag 1 in sigma1-cyclische ordening
    f0_c  = f0 - f0.mean()
    f0s_c = f0_shifted - f0_shifted.mean()   # == f0_shifted - f0.mean() (sigma1 is bijectie)

    cov_01 = float(np.mean(f0_c * f0s_c))
    var_f0 = float(np.var(f0))
    rho1   = cov_01 / var_f0 if var_f0 > 0 else 0.0

    # Verify: Cov[f0(s), f1(s)] / Var[f0] should ≈ rho1
    # (want v1(s) = (A/rho)*v0(sigma1(s)), dus f1 = f0 circ sigma1 + const)
    f1_c = f1 - f1.mean()
    cov_f0f1 = float(np.mean(f0_c * f1_c))
    rho_cross = cov_f0f1 / var_f0 if var_f0 > 0 else 0.0

    # Mixing gain: Var[f0(sigma1(s)) - f0(s)] / Var[f0] = 2*(1 - rho1)
    diff = f0_shifted - f0
    mixing_gain = float(np.var(diff)) / var_f0 if var_f0 > 0 else 0.0

    # r=1 check: ve1 / ve0 gives the local ratio; compare with rho1
    # Naieve verwachting: ve1 ~ ve0 * rho1^2? (nee, via X1 = f0(sigma1) - lm)
    # Betere vergelijking: ve1 / var_f0 als functie van rho1

    print(f"  var_end       = {var_end:.8f}")
    print(f"  ve0, ve1, ve2 = {ve0:.6f}, {ve1:.6f}, {ve2:.6f}")
    print(f"  (ve0+ve1+ve2)/3 = {(ve0+ve1+ve2)/3:.8f}  (check ~= var_end)")
    print(f"  rho1 = Cov[f0(s), f0(sigma1(s))] / Var[f0] = {rho1:.6f}")
    print(f"  rho_cross = Cov[f0(s), f1(s)] / Var[f0]    = {rho_cross:.6f}  (check: ~= rho1)")
    print(f"  mixing_gain = Var[f0(sigma1)-f0] / Var[f0]  = {mixing_gain:.6f}  (check: ~= 2*(1-rho1)={2*(1-rho1):.6f})")
    print()
    sys.stdout.flush()

    results[k] = {
        'var_end': var_end, 've0': ve0, 've1': ve1, 've2': ve2,
        'rho1': rho1, 'rho_cross': rho_cross, 'mixing_gain': mixing_gain
    }

# Samenvatting
print("=" * 70)
print("SAMENVATTING")
print()
print(f"{'k':>3}  {'var_end':>10}  {'d_k-1':>8}  {'rho1':>8}  {'mix_gain':>9}  {'ve1/ve0':>8}")
print("-" * 65)
prev_ve = None
for k in K_RANGE:
    r = results[k]
    dk_prev = r['var_end'] / prev_ve if prev_ve else float('nan')
    ve1_ve0 = r['ve1'] / r['ve0'] if r['ve0'] > 0 else float('nan')
    print(f"{k:3d}  {r['var_end']:10.8f}  {dk_prev:8.6f}  {r['rho1']:8.6f}  {r['mixing_gain']:9.6f}  {ve1_ve0:8.6f}")
    prev_ve = r['var_end']

print()
print("Theoretische verwachting:")
print("  rho1 < 1  =>  sigma1-shuffle introduceert diversiteit => bijdrage aan d_k < 1")
print("  mixing_gain = 2*(1-rho1): de relatieve variantie van de 'shift'")
print("  ve1 < ve0: door sigma1-mixing krimpt de r=1 bijdrage t.o.v. r=0")
print()
print("done")
