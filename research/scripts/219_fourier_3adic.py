"""
219_fourier_3adic.py
====================
Methode 10: Fourier-analyse op Z/3^k Z van de Perron-eigenvector.

De Perron-vector v^(k) leeft op Z/3^{k-1}Z.
Haar Fourier-transformatie op deze groep is:
  v̂(n) = (1/N) Σ_{j=0}^{N-1} v(j) * omega^{-jn}
waarbij omega = exp(2*pi*i/N) een primitieve N-de eenheidswortel is.

De 3-adische Hoelder-conditie (onze karakterisering van Conjecture G)
is equivalent met:
  |v̂(n)| <= C * |n|_3^{alpha}   voor zekere alpha > 0

waarbij |n|_3 = 3^{-v_3(n)} de 3-adische norm is van n.

In termen van de DFT-indices: n met hoge 3-adische valuatie v_3(n)=l
correspondeert met "lage 3-adische frequenties" (structuur op groot schaal).
n met lage v_3(n) = hoge 3-adische frequenties (fijne structuur).

De Hoelder-conditie zegt: de DFT-coëfficiënten vallen AF als functie
van de 3-adische norm |n|_3.

We meten:
  (A) De DFT van v^(k) (grootte en fase)
  (B) Afvalgedrag van |v̂(n)| als functie van |n|_3
  (C) Helling van log|v̂(n)| vs log|n|_3 (= Hoelder exponent alpha)
  (D) Of alpha > 0 (stijgend bewijs voor G)

Als alpha_k -> alpha_inf > 0 als k -> inf: direct bewijs voor G.
"""
import numpy as np
from math import log2, log

ALPHA_LOG = log2(3.0)
LAM   = 1.70
A     = LAM ** -2.0
B1    = LAM ** (ALPHA_LOG - 2.0)
B3    = LAM ** (ALPHA_LOG - 1.0)


def perron(k, n_iter=300):
    N  = 3 ** (k - 1)
    i  = np.arange(N, dtype=np.int64)
    T4 = (4 * i + 2) % N
    s, r = np.divmod(i, 3)
    Nl = N // 3
    m0, m2 = (r == 0), (r == 2)
    R1 = (4 * s) % Nl
    R3 = (2 * s + 1) % Nl
    v = np.ones(N, dtype=np.float64)
    for _ in range(n_iter):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        w  = A * v[T4]
        w[m2] += B3 * cb[R3[m2]]
        w[m0] += B1 * cb[R1[m0]]
        v = w / w.max()
    v /= v.mean()
    return v


def v3_val(n):
    """3-adische valuatie van n (0 als n=0)."""
    if n == 0:
        return 100
    cnt = 0
    while n % 3 == 0:
        n //= 3
        cnt += 1
    return cnt


def analyse_fourier(k):
    v = perron(k)
    N = v.size

    # DFT
    vhat = np.fft.fft(v) / N

    # Voor elk n: |v̂(n)| en |n|_3 = 3^{-v3(n)}
    ns   = np.arange(N, dtype=np.int64)
    mags = np.abs(vhat)

    # v_3(n) voor elk n (0 heeft valuatie oneindig -> sluit 0 uit)
    v3s  = np.array([v3_val(int(n)) for n in ns[1:]], dtype=np.float64)
    mag1 = mags[1:]  # DFT buiten de DC-component

    # 3-adische norm |n|_3 = 3^{-v3(n)}
    norm3 = 3.0 ** (-v3s)

    # Groepeer per 3-adische valuatielaag
    max_v3 = int(v3s.max())
    print(f"  k={k}  N={N}", flush=True)
    print(f"  {'v3':>4s}  {'|n|_3':>8s}  {'#n':>6s}  "
          f"{'mean|vhat|':>10s}  {'max|vhat|':>10s}", flush=True)

    layer_means = []
    layer_norms = []
    for l in range(0, min(max_v3 + 1, k)):
        mask = (v3s == l)
        if mask.sum() == 0:
            continue
        mean_mag = float(mag1[mask].mean())
        max_mag  = float(mag1[mask].max())
        norm_l   = 3.0 ** (-l)
        print(f"  {l:>4d}  {norm_l:>8.5f}  {mask.sum():>6d}  "
              f"{mean_mag:>10.6f}  {max_mag:>10.6f}", flush=True)
        layer_means.append(mean_mag)
        layer_norms.append(norm_l)

    # Fit log|v̂| vs log|n|_3 (Hoelder exponent)
    layer_means = np.array(layer_means)
    layer_norms = np.array(layer_norms)
    ok = layer_means > 1e-12
    if ok.sum() >= 3:
        log_norms = np.log(layer_norms[ok])
        log_mags  = np.log(layer_means[ok])
        slope, intercept = np.polyfit(log_norms, log_mags, 1)
        print(f"  Hoelder exponent alpha = {slope:.4f}  "
              f"(log|v̂| = {slope:.4f} * log|n|_3 + {intercept:.4f})", flush=True)
        print(f"  Interpretatie: alpha > 0 -> Holder-regulier -> G bevestigd",
              flush=True)
        return slope
    return float('nan')


print(f"Methode 10: Fourier-analyse op Z/3^kZ  (lam={LAM})", flush=True)
print("=" * 65, flush=True)

alphas = {}
for k in (10, 12, 14, 15, 16):
    print(flush=True)
    alpha = analyse_fourier(k)
    alphas[k] = alpha

print(f"\nHoelder exponent samenvatting:", flush=True)
ks = sorted(alphas.keys())
for k in ks:
    print(f"  k={k:2d}:  alpha = {alphas[k]:.5f}", flush=True)

# Extrapolatie: convergeert alpha_k naar een positieve limiet?
if len(ks) >= 3:
    ks_arr = np.array(ks, dtype=float)
    al_arr = np.array([alphas[k] for k in ks])
    ok = np.isfinite(al_arr)
    if ok.sum() >= 3:
        slope_alpha, intercept_alpha = np.polyfit(ks_arr[ok], al_arr[ok], 1)
        alpha_inf = intercept_alpha + slope_alpha * 100  # extrapolatie k=100
        print(f"\n  Lineaire extrapolatie: alpha(k) = {slope_alpha:.5f}*k + "
              f"{intercept_alpha:.5f}", flush=True)
        print(f"  alpha_inf (k->inf) ~ {alpha_inf:.4f}", flush=True)
        print(f"  alpha_inf > 0: {'JA -> G waarschijnlijk' if alpha_inf > 0 else 'NEE'}", flush=True)

print("\ndone", flush=True)
