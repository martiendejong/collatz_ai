"""
227_fourier_hoelder_k17.py
===========================
Fourier-Hoelder-exponent van de K-L Perron-eigenvector voor k=17,18.

Verlenging van Obs 424 (k=10,12,14,15,16).
Doel: aanscherpen van alpha_inf >= 0.646 (Conjecture G numeriek).

Methode: DFT van v^(k) op Z/3^{k-1}Z. Bereken |vhat(n)| voor alle n.
Groepeer per 3-adische valuatie v3(n): v3(n) = max{j: 3^j | n}.
Log-lineaire regressie van log|vhat| vs log|n|_3 geeft Hoelder-exponent alpha.

Eerdere resultaten:
  k=10: |alpha| = 0.706
  k=12: |alpha| = 0.687
  k=14: |alpha| = 0.675
  k=15: |alpha| = 0.670
  k=16: |alpha| = 0.666

Extrapolatie: decrements krimpen geometrisch (ratio ~0.8).
Resterende serie vanaf k=16: 0.004/(1-0.8) = 0.020.
Huidige bound: alpha_inf >= 0.646.

Doel van deze script: k=17 en eventueel k=18 berekenen om de extrapolatie
te verfijnen.

NOOT: k=17 geeft N = 3^16 = 43046721 knopen. DFT hiervan is O(N log N)
en vereist ~1 GB geheugen. Dit kan ZWAAR zijn.
Voor k=18: N = 3^17 = 129140163. Mogelijk te groot.

Aanpak: bereken alleen de DFT voor de subgroep van knopen met v3-index
in {0, 1, ..., k-2} (= alle knopen behalve n=0).
Bereken |vhat(n)| via directe sommatie voor een STEEKPROEF van n-waarden.
"""
import numpy as np
from math import log2, log

ALPHA_MATH = log2(3.0)
LAM   = 1.70
A  = LAM ** -2.0
B1 = LAM ** (ALPHA_MATH - 2.0)
B3 = LAM ** (ALPHA_MATH - 1.0)


def three_adic_valuation(n, N):
    """v3(n) = max{j: 3^j | n}, met wraparound mod N."""
    n = int(n) % int(N)
    if n == 0:
        return 999  # de nul-mode
    v = 0
    while n % 3 == 0:
        n //= 3; v += 1
    return v


def compute_alpha_k(k, n_iter=400, n_samples=2000):
    """
    Bereken Hoelder-exponent alpha voor grootte k.
    Gebruik steekproef van n-waarden om te vermijden dat we de volledige DFT nodig hebben.
    """
    N  = 3 ** (k - 1)
    i  = np.arange(N, dtype=np.int64)
    T4 = (4 * i + 2) % N
    s, r = np.divmod(i, 3)
    Nl = N // 3
    m0, m2 = (r == 0), (r == 2)
    R1 = (4 * s) % Nl
    R3 = (2 * s + 1) % Nl

    # Perron-eigenvector (gemiddelde = 1)
    v = np.ones(N, dtype=np.float64)
    for _ in range(n_iter):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        w  = A * v[T4]
        w[m2] += B3 * cb[R3[m2]]
        w[m0] += B1 * cb[R1[m0]]
        v = w / w.max()
    v /= v.mean()

    print(f"  k={k}: N={N}, eigenvector berekend.")

    # Bereken DFT via numpy FFT als het past, anders via steekproef.
    # Voor k<=14: N<=1594323, FFT is OK.
    # Voor k=15: N=4782969 (~5M), FFT gaat ~4 sec.
    # Voor k=16: N=14348907 (~14M), FFT gaat ~15 sec.
    # Voor k=17: N=43046721 (~43M), FFT gaat ~50 sec en ~350MB RAM.
    # Voor k=18: N=129140163 (~129M), te groot.

    if N <= 50_000_000:
        # Volledige FFT
        vhat = np.fft.fft(v)
        abs_vhat = np.abs(vhat)

        # Groepeer per 3-adische valuatie
        # v3(n) voor n in {1,...,N-1}
        # Gemiddeld |vhat(n)| per v3-niveau
        v3_levels = {}
        for n in range(1, N):
            v3 = 0
            m = n
            while m % 3 == 0:
                m //= 3; v3 += 1
            if v3 not in v3_levels:
                v3_levels[v3] = []
            v3_levels[v3].append(abs_vhat[n])

        # Gemiddeld log|vhat| per niveau
        pts = []
        for v3, vals in sorted(v3_levels.items()):
            mean_log_vhat = float(np.mean(np.log(np.array(vals) + 1e-300)))
            # |n|_3 = 3^{-v3} (de 3-adische norm)
            log_norm = -v3 * log(3.0)   # log(3^{-v3}) = -v3*log(3)
            pts.append((log_norm, mean_log_vhat))

        pts = sorted(pts)
        log_norms = np.array([p[0] for p in pts])
        log_vhats = np.array([p[1] for p in pts])

        # Lineaire regressie (Hoelder-exponent = helling)
        ok = np.isfinite(log_vhats) & (log_norms > -15)
        if ok.sum() >= 3:
            slope, intercept = np.polyfit(log_norms[ok], log_vhats[ok], 1)
            alpha = slope   # alpha = d(log|vhat|)/d(log|n|_3)
        else:
            slope = float('nan')
            alpha = float('nan')

        print(f"  k={k}: alpha = {alpha:.4f} (via FFT, {len(pts)} niveaus)")
        # Toon ook de ruwe data
        print(f"  Niveaus: {[(round(p[0],2), round(p[1],3)) for p in pts[:8]]}")
        return alpha

    else:
        # Steekproef: kies n met v3(n) = 0,1,2,...,k-3 elk vertegenwoordigd
        print(f"  k={k}: N te groot voor FFT ({N}), gebruik steekproef.")
        # Kies representatieve n per v3-niveau
        means = {}
        for target_v3 in range(k - 2):
            # Vind ~100 n-waarden met v3(n) = target_v3
            # n = 3^target_v3 * m met gcd(m,3)=1
            base = 3 ** target_v3
            candidates = [base * m for m in range(1, N // base, (N // base) // 100 + 1)
                          if m % 3 != 0 and base * m < N][:100]
            if not candidates:
                continue
            # Bereken DFT bij deze n
            vhat_vals = []
            for n_val in candidates:
                phases = np.exp(-2j * np.pi * n_val * np.arange(N) / N)
                vhat_n = float(np.abs(np.dot(v, phases)))
                vhat_vals.append(vhat_n)
            mean_log = float(np.mean(np.log(np.array(vhat_vals) + 1e-300)))
            means[target_v3] = mean_log

        pts = []
        for v3, mean_log in sorted(means.items()):
            log_norm = -v3 * log(3.0)
            pts.append((log_norm, mean_log))

        log_norms = np.array([p[0] for p in pts])
        log_vhats = np.array([p[1] for p in pts])
        ok = np.isfinite(log_vhats)
        if ok.sum() >= 3:
            slope, _ = np.polyfit(log_norms[ok], log_vhats[ok], 1)
            alpha = slope
        else:
            alpha = float('nan')

        print(f"  k={k}: alpha = {alpha:.4f} (via steekproef)")
        return alpha


print(f"227: Fourier-Hoelder-exponent k=14..17  (lam={LAM})")
print("=" * 60)
print()
print("Eerdere resultaten (Obs 424):")
print("  k=10: 0.706, k=12: 0.687, k=14: 0.675, k=15: 0.670, k=16: 0.666")
print()

results = {}
for k in range(14, 18):
    alpha = compute_alpha_k(k)
    results[k] = alpha
    print()

print("=== SAMENVATTING ===")
ks = sorted(results.keys())
alphas = [results[k] for k in ks]
print(f"k     alpha")
for k, a in zip(ks, alphas):
    print(f"  {k:2d}  {a:.4f}")

print()
# Decrements en extrapolatie
decrements = np.diff(alphas)
print(f"Decrements: {[f'{d:.4f}' for d in decrements]}")
if len(decrements) >= 2:
    ratios = decrements[1:] / np.abs(decrements[:-1] + 1e-12)
    print(f"Ratios: {[f'{r:.3f}' for r in ratios]}")
    last_dec = abs(decrements[-1])
    last_ratio = abs(ratios[-1])
    if last_ratio < 1:
        tail_sum = last_dec / (1 - last_ratio)
        alpha_inf_lower = alphas[-1] - tail_sum
        print(f"Geextrapoleerde alpha_inf >= {alpha_inf_lower:.4f}")
    else:
        print("Ratio >= 1: extrapolatie niet betrouwbaar")

print()
print("done")
