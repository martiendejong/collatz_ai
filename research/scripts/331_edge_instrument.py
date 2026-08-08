"""
Test: is the ~1.10 cross-instrument factor the edge layer?
Compute per k (lam=1.70; k=12,13,14): the full inc-spectrum, then
  (a) across-k ratio of the TOP-layer increment inc_{k-2}(k) — the d-analog;
  (b) across-k ratio of DEEP-TAIL sums T_e(k) = sum of last e increments, e=1..4;
  (c) the within-k bulk plateau (edge-excluded) — the r-instrument.
If top-layer/tail ratios across k ~ 0.76 (=1.10*r) while bulk plateau ~0.69,
the edge-layer explanation of d/r ~ 1.10 is confirmed.
"""
import numpy as np
from math import log2
CACHE = "E:/projects/collatz/research/cache"

def spectrum(lam, k):
    v = np.load(f"{CACHE}/v_lam{lam:.2f}_k{k}.npy")
    N = v.size
    F = np.log2(v); F -= F.mean()
    incs = []; prev = None
    for p in range(k-1):
        M = 3**(p+1)
        cm = F.reshape(N//M, M).mean(axis=0)
        m = cm if prev is None else cm - prev[np.arange(M) % (M//3)]
        incs.append(float((m**2).mean())); prev = cm
    return np.array(incs)

lam = 1.70
S = {k: spectrum(lam, k) for k in [12, 13, 14]}
print("inc-spectra tails (last 5) per k:")
for k in [12, 13, 14]:
    print(f"  k={k}:", " ".join(f"{x:.3e}" for x in S[k][-5:]))

print("\n(a) top-layer increment across k: inc_{k-2}(k)")
tops = {k: S[k][-1] for k in [12, 13, 14]}
print(f"  values: {tops[12]:.3e} {tops[13]:.3e} {tops[14]:.3e}")
print(f"  ratios: {tops[13]/tops[12]:.4f} {tops[14]/tops[13]:.4f}")

print("\n(b) deep-tail sums T_e across k (e = laatste e lagen):")
for e in [1, 2, 3, 4]:
    T = {k: S[k][-e:].sum() for k in [12, 13, 14]}
    print(f"  e={e}: ratios {T[13]/T[12]:.4f} {T[14]/T[13]:.4f}")

print("\n(c) binnen-k bulk plateau (randlagen uitgesloten, p=4..k-5):")
for k in [12, 13, 14]:
    rats = S[k][1:]/S[k][:-1]
    bulk = rats[4:k-5]
    print(f"  k={k}: bulk median {np.median(bulk):.4f}  (edge ratios: {' '.join(f'{x:.3f}' for x in rats[-3:])})")
print("\npaper d-plateau at lam=1.70: 0.756-0.766 (k=13-14); r-bulk ~0.69 => d/r ~ 1.10")
