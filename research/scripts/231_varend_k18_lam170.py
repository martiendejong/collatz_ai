"""
231_varend_k18_lam170.py
========================
Bereken var_end(18, lambda=1.70) via directe power-iteratie.
N=3^17=129,140,163. 200 iteraties. Memory: ~4 GB float64.

Geeft:
  d_17(lam=1.70) = var_end(18)/var_end(17) = var_end(18)/0.00065480
  d_18(lam=1.70) = var_end(19)/var_end(18) = 0.00038900/var_end(18)

Bekende waarden:
  var_end(17, 1.70) = 0.00065480  (Script 229)
  var_end(19, 1.70) = 0.00038900  (Script 200b)
  var_end(20, 1.70) = 0.00030161  (Script 229b) -> d_19=0.7753
"""
import sys
import numpy as np
from math import log2

ALPHA = log2(3.0)
LAM   = 1.70
A  = LAM ** -2.0
B1 = LAM ** (ALPHA - 2.0)
B3 = LAM ** (ALPHA - 1.0)
N_ITER = 200

VAREND_17 = 0.00065480   # Script 229
VAREND_19 = 0.00038900   # Script 200b

print(f"231: var_end(18, lam={LAM}) directe berekening")
print(f"     A={A:.6f}  B1={B1:.6f}  B3={B3:.6f}")
print(f"     N_iter={N_ITER}")
print("=" * 60)
sys.stdout.flush()

k = 18
N  = 3 ** (k - 1)
Nl = N // 3
print(f"  k={k}: N={N:,}  Nl={Nl:,}")
print(f"  Memory estimate: {N * 8 / 1e9:.2f} GB (float64)")
sys.stdout.flush()

i  = np.arange(N, dtype=np.int64)
T4 = (4 * i + 2) % N
s, r = np.divmod(i, 3)
m0, m2 = (r == 0), (r == 2)
R1 = (4 * s) % Nl
R3 = (2 * s + 1) % Nl

print("  Indexen berekend. Start power-iteratie ...")
sys.stdout.flush()

v = np.ones(N, dtype=np.float64)

for it in range(N_ITER):
    cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
    w  = A * v[T4]
    w[m2] += B3 * cb[R3[m2]]
    w[m0] += B1 * cb[R1[m0]]
    w_max = w.max()
    v = w / w_max

    if (it + 1) % 50 == 0 or it == N_ITER - 1:
        # rho check
        cb2 = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        w2 = A * v[T4]
        w2[m2] += B3 * cb2[R3[m2]]
        w2[m0] += B1 * cb2[R1[m0]]
        rho = float(w2.mean() / v.mean())
        print(f"  iter={it+1:3d}: rho={rho:.6f}")
        sys.stdout.flush()

print()

# var_end berekening
print("  var_end berekening ...")
T = np.stack([v[:Nl], v[Nl:2*Nl], v[2*Nl:]])   # (3, Nl)
lmean = T.mean(axis=0)
X = np.log2(T) - np.log2(lmean)[None, :]
ve18 = float(np.var(X))
print(f"  var_end(18, lam=1.70) = {ve18:.8f}")
sys.stdout.flush()

# d_k afgeleid
d17 = ve18 / VAREND_17
d18 = VAREND_19 / ve18
print()
print("=== RESULTATEN ===")
print(f"  var_end(17) = {VAREND_17:.8f}  (Script 229)")
print(f"  var_end(18) = {ve18:.8f}  (dit script)")
print(f"  var_end(19) = {VAREND_19:.8f}  (Script 200b)")
print(f"  d_17 = var_end(18)/var_end(17) = {d17:.6f}")
print(f"  d_18 = var_end(19)/var_end(18) = {d18:.6f}")
print()
print("Volledige d_k-reeks (lambda=1.70) ZONDER GAP:")
dk_all = {
    13: 0.756036, 14: 0.753544, 15: 0.759039,
    16: 0.766188, 17: d17, 18: d18, 19: 0.775346
}
for kk, dk in sorted(dk_all.items()):
    print(f"  d_{kk} = {dk:.6f}")

dk_vals = list(dk_all.values())
print(f"\nGemiddeld d_k (k=13..19): {np.mean(dk_vals):.5f}")
print(f"Alle d_k < 1: {'JA' if all(d < 1 for d in dk_vals) else 'NEE'}")
print(f"Trend increment d_17->d_18->d_19: {d17-dk_all[16]:.5f}, {d18-d17:.5f}, {0.775346-d18:.5f}")
print()
print("done")
