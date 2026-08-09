"""
The linearization test: is the within-triple deviation shape asymptotically Gaussian?
For sum-zero Gaussian triples the shape constant kappa = E[mean-min]/E[sigma_w] is
universal (computable by simulation). Measure at the fixed point across lambda and
depth-slices; compare moments (skew, kurtosis of normalized deviations) with Gaussian.
If Gaussian: the min becomes 'linear + bilinear with known kappa' — the linearization.
"""
import numpy as np
from math import log2
CACHE = "E:/projects/collatz/research/cache"

# Gaussian reference for sum-zero triples
rng = np.random.default_rng(3)
G = rng.standard_normal((2_000_000, 3))
G = G - G.mean(axis=1, keepdims=True)
sw = G.std(axis=1)
kappa_gauss = float((G.mean(axis=1) - G.min(axis=1)).mean() / sw.mean())
gap_g = -G.min(axis=1)
print(f"Gauss-referentie: kappa = E[gap]/E[sigma_w] = {kappa_gauss:.5f}")
d_g = (G / sw[:, None]).ravel()
print(f"  genormaliseerde deviaties: skew {float(((d_g-d_g.mean())**3).mean()):+.4f} "
      f"kurt-3 {float(((d_g-d_g.mean())**4).mean())-3:+.4f}\n")

for lam in [1.05, 1.70, 2.00]:
    k = 13
    v = np.load(f"{CACHE}/v_lam{lam:.2f}_k{k}.npy")
    N = v.size; Nl = N//3
    trip = np.stack([v[:Nl], v[Nl:2*Nl], v[2*Nl:]], axis=1)
    dev = trip - trip.mean(axis=1, keepdims=True)
    sw = dev.std(axis=1)
    ok = sw > 0
    kappa = float((-dev.min(axis=1))[ok].mean() / sw[ok].mean())
    nd = (dev[ok] / sw[ok, None]).ravel()
    skew = float(((nd - nd.mean())**3).mean())
    kurt = float(((nd - nd.mean())**4).mean()) - 3
    print(f"lam={lam}: kappa = {kappa:.5f} (Gauss {kappa_gauss:.5f}, ratio {kappa/kappa_gauss:.4f}) "
          f"| skew {skew:+.4f} kurt-3 {kurt:+.4f}")
