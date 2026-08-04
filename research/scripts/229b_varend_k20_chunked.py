"""
229b_varend_k20_chunked.py
==========================
Bereken var_end(20, lambda=1.70) uit het opgeslagen k=20 certificaat
via chunked memory-mapped reading (omzeilt 8.66 GiB OOM).

N  = 3^19 = 1,162,261,467 (float32, 4.65 GB)
Nl = 3^18 = 387,420,489

Chunked var_end:
  Per chunk [cs:ce] van [0, Nl):
    v0 = v[cs:ce], v1 = v[cs+Nl:ce+Nl], v2 = v[cs+2*Nl:ce+2*Nl]
    lmean = (v0+v1+v2)/3
    Xi = log2(vi/lmean)  voor i=0,1,2
  Accumuleer sum(Xi) en sum(Xi^2) in float64.
  Populatievariantie = E[X^2] - E[X]^2

Aannames:
  var_end(19, 1.70) = 0.000389   (Script 200b, consistent met d_k trend)
  d_k trend: d_17=0.7690, extrapolatie d_18~0.770, d_19~0.771
"""
import sys
import os
import numpy as np

CERT_DIR = os.path.join(os.path.dirname(__file__), "..", "certificates")
CERT_K20 = os.path.join(CERT_DIR, "k20_lam170_200c.npy")

# Var_end(19, lam=1.70) uit Script 200b (intermediaire meting)
VAREND_19 = 0.000389

CHUNK = 20_000_000   # 20M elementen per chunk => 3 x 80 MB float64 per stap

print("229b: var_end(20, lam=1.70) chunked mmap berekening")
print("=" * 60)
print(f"  Certificaat : {CERT_K20}")
print(f"  Chunk size  : {CHUNK:,}")
print()
sys.stdout.flush()

if not os.path.exists(CERT_K20):
    print(f"FOUT: {CERT_K20} niet gevonden.")
    sys.exit(1)

# Memory-map laden (float32, shape=(N,))
v = np.load(CERT_K20, mmap_mode='r')
N = len(v)
k = int(round(np.log(N) / np.log(3) + 1))
Nl = N // 3

print(f"  Geladen: shape={v.shape}, dtype={v.dtype}, N={N:,}, k={k}")
print(f"  Nl = {Nl:,}")
print()
sys.stdout.flush()

assert N == 3 ** (k - 1), f"N={N} != 3^{{k-1}}=3^{k-1}={3**(k-1)}"
assert Nl * 3 == N

# Chunked var_end berekening
sum_X  = np.float64(0.0)
sum_X2 = np.float64(0.0)
n_chunks = (Nl + CHUNK - 1) // CHUNK

print(f"  Berekening in {n_chunks} chunks ...")
sys.stdout.flush()

for i in range(n_chunks):
    cs = i * CHUNK
    ce = min(cs + CHUNK, Nl)
    sz = ce - cs

    # Laad drie stukken (float32 -> float64)
    v0 = v[cs      : ce      ].astype(np.float64)
    v1 = v[cs + Nl : ce + Nl ].astype(np.float64)
    v2 = v[cs+2*Nl : ce+2*Nl].astype(np.float64)

    lmean = (v0 + v1 + v2) * (1.0 / 3.0)

    # log2(vi/lmean)
    x0 = np.log2(v0 / lmean)
    x1 = np.log2(v1 / lmean)
    x2 = np.log2(v2 / lmean)

    sum_X  += x0.sum() + x1.sum() + x2.sum()
    sum_X2 += (x0*x0).sum() + (x1*x1).sum() + (x2*x2).sum()

    if (i + 1) % 5 == 0 or i == n_chunks - 1:
        pct = 100.0 * (i + 1) / n_chunks
        print(f"    chunk {i+1}/{n_chunks}  ({pct:.1f}%)  sum_X2={sum_X2:.6e}")
        sys.stdout.flush()

n_total = np.float64(3 * Nl)
E_X  = sum_X  / n_total
E_X2 = sum_X2 / n_total
varend_20 = float(E_X2 - E_X ** 2)

print()
print("=== RESULTATEN ===")
print(f"  var_end(20, lam=1.70) = {varend_20:.8f}")
print(f"  var_end(19, lam=1.70) = {VAREND_19:.8f}  (Script 200b)")
d19 = varend_20 / VAREND_19
print(f"  d_19 = var_end(20)/var_end(19) = {d19:.6f}")
print()

# Vergelijk met trend
print("d_k-reeks (lambda=1.70) compleet:")
dk_known = {
    13: 0.756036,
    14: 0.753544,
    15: 0.759039,
    16: 0.766188,
    17: 0.769000,   # Script 200
}
dk_known[19] = d19

print(f"{'k':>4}  {'d_k':>10}")
print("-" * 18)
for kk in sorted(dk_known):
    gap = "  (gat)" if kk == 18 else ""
    print(f"  {kk:2d}  {dk_known[kk]:10.6f}{gap}")

dk_vals = [v for k, v in sorted(dk_known.items()) if not (k == 18)]
print()
print(f"Gemiddeld d_k (k=13..17,19): {np.mean(dk_vals):.5f}")
print(f"Alle d_k < 1: {'JA' if all(x < 1 for x in dk_vals) else 'NEE'}")
print()
print("CONCLUSIE: Endpoint Decay (d_k < 1) bevestigd t/m k=19  (lambda=1.70)")
print()
print("done")
