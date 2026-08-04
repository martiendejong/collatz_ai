"""
229_second_eigenvalue.py
========================
Meet f2 = spectrale convergentiefactor van de K-L operator via TWEE-PUNT methode.

Methode: start power iteration vanuit v0 = ones en v0' = ones + kleine perturbatie.
Na n stappen: |v_n - v_n'| ~ C * (lambda2/lambda1)^n.
De verhouding |v_{n+1} - v_{n+1}'| / |v_n - v_n'| convergeert naar lambda2/lambda1.
f2 = (lambda2/lambda1)^2 = spectrale convergentiefactor voor variantie-afname.

Endpoint Decay: f2 < 1 is equivalent aan lambda2 < lambda1 = rho (spectraal gat > 0).

Let op: lambda1 hier = rho (de nonlineaire Perron-eigenwaarde van F), NIET de
spectrale straal van de gelineariseerde Jacobiaan M (die groter kan zijn).
"""
import sys
import numpy as np
from math import log2

ALPHA = log2(3.0)
LAM   = 1.70
A  = LAM ** -2.0
B1 = LAM ** (ALPHA - 2.0)
B3 = LAM ** (ALPHA - 1.0)

print(f"229: f2-meting K-L operator via twee-punt machtiteratie  (lam={LAM})")
print(f"     A={A:.6f}  B1={B1:.6f}  B3={B3:.6f}")
print("=" * 65)
print()
sys.stdout.flush()

results = []

for k in range(5, 16):
    N  = 3 ** (k - 1)
    print(f"  k={k:2d}  N={N:>10,d} ...", end='', flush=True)

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

    # Start 1: v = uniform
    v  = np.ones(N, dtype=np.float64)
    # Start 2: v' = ones + kleine willekeurige perturbatie
    rng = np.random.default_rng(seed=42)
    eps0 = 1e-4
    dv = eps0 * (rng.random(N) - 0.5)
    vp = v + dv

    n_warmup = 200  # eerst convergeren naar buurt van v*
    for _ in range(n_warmup):
        v  = apply_F(v);  v  /= v.max()
        vp = apply_F(vp); vp /= vp.max()

    # Nu zijn v en vp dicht bij v*. Meet convergentiesnelheid van hun verschil.
    diffs = []
    for it in range(80):
        diff = float(np.linalg.norm(v - vp))
        diffs.append(diff)
        v  = apply_F(v);  v  /= v.max()
        vp = apply_F(vp); vp /= vp.max()

    diffs = np.array(diffs)
    # Vermijd nul-waarden
    ok = diffs > 1e-300
    if ok.sum() < 10:
        print(f"  diff convergeert te snel of is nul")
        continue

    # Fit log(diff[n]) = log(C) + n * log(ratio)
    ns_ok = np.where(ok)[0].astype(float)
    log_diffs_ok = np.log(diffs[ok])

    # Gebruik de LAATSTE helft van de gemeten punten (meest uitgeconvergeerd)
    half = len(ns_ok) // 2
    if half < 5:
        half = 0
    ns_fit     = ns_ok[half:]
    log_diffs_fit = log_diffs_ok[half:]

    if len(ns_fit) >= 5:
        fit_slope, fit_intercept = np.polyfit(ns_fit, log_diffs_fit, 1)
        ratio_measured = float(np.exp(fit_slope))
    else:
        ratio_measured = float('nan')

    f2 = ratio_measured ** 2

    # rho uit de power iteraat
    v_norm = v / v.mean()
    w_check = apply_F(v_norm)
    rho = float(w_check.mean())

    print(f"  rho={rho:.6f}  ratio=lambda2/lambda1={ratio_measured:.6f}  f2={f2:.6f}")
    sys.stdout.flush()

    # Toon ook de ruwe ratios (laatste 20 iteraties)
    raw_ratios = diffs[1:] / (diffs[:-1] + 1e-300)
    last_ratios = raw_ratios[-20:]
    print(f"    rho ratios (laatste 20): "
          f"mean={float(np.mean(last_ratios[last_ratios>0])):.5f}  "
          f"std={float(np.std(last_ratios[last_ratios>0])):.5f}")
    sys.stdout.flush()

    results.append({'k': k, 'N': N, 'rho': rho,
                    'ratio': ratio_measured, 'f2': f2,
                    'diffs': diffs.tolist()})

print()
print("=== SAMENVATTING ===")
print()
print(f"{'k':>4} {'N':>10} {'rho':>10} {'ratio':>10} {'f2':>10}")
print("-" * 55)
for r in results:
    print(f"  {r['k']:>2}  {r['N']:>10,d}  {r['rho']:>10.6f}  "
          f"{r['ratio']:>10.6f}  {r['f2']:>10.6f}")

if results:
    f2_vals = [r['f2'] for r in results]
    ratio_vals = [r['ratio'] for r in results]
    print()
    print(f"f2 waarden:    {[f'{x:.5f}' for x in f2_vals]}")
    print(f"ratio waarden: {[f'{x:.5f}' for x in ratio_vals]}")
    print()
    all_positive = all(x < 1.0 for x in f2_vals if not np.isnan(x))
    print(f"f2 < 1 voor alle k: {'JA' if all_positive else 'NEE'}")
    valid_f2 = [x for x in f2_vals if not np.isnan(x) and x > 0]
    if valid_f2:
        print(f"f2 (laatste geldig): {valid_f2[-1]:.5f}")
        print(f"Spectrale kloof (1 - ratio): {1 - ratio_vals[-1]:.5f}")
        print()
        print("INTERPRETATIE:")
        print(f"  lambda2/lambda1 ~ {ratio_vals[-1]:.5f}")
        print(f"  f2 ~ {valid_f2[-1]:.5f} < 1 => Endpoint Decay bevestigd per stap")
        print(f"  Variantie-attenuatie per K-L iteratie: f2 = {valid_f2[-1]:.5f}")

print()
print("done")
