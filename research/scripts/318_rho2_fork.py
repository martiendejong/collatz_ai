"""
Fork discriminator: 1 - rho(lambda=2, k) for k=5..16.
DENSITY model: rho(2,inf)=1 -> pure geometric decay to 0, rate = endpoint contraction c.
CEILING model: rho(2,inf)<1 -> decay levels off at a positive constant.
Fit both: y_k = C*r^k  vs  y_k = y_inf + C*r^k, compare residuals.
"""
import numpy as np
from math import log2

ALPHA = log2(3.0)
lam = 2.0
A = lam**-2; B1 = lam**(ALPHA-2); B3 = lam**(ALPHA-1)
ITERS = {5:4000, 6:4000, 7:3000, 8:3000, 9:2500, 10:2000, 11:1500, 12:1200, 13:900, 14:600, 15:400, 16:300}

rhos = {}
for k in range(5, 17):
    N = 3**(k-1); Nl = N//3
    i = np.arange(N, dtype=np.int64)
    T4 = (4*i+2) % N
    s_arr, r_arr = np.divmod(i, 3)
    m0 = r_arr==0; m2 = r_arr==2
    R1 = (4*s_arr) % Nl; R3 = (2*s_arr+1) % Nl
    del i, s_arr, r_arr
    v = np.ones(N)
    rho = 1.0
    for _ in range(ITERS[k]):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        w = A*v[T4]
        w[m2] += B3*cb[R3[m2]]
        w[m0] += B1*cb[R1[m0]]
        rho = float(w.max()); w /= rho; v = w
    # CW check
    cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
    w = A*v[T4]; w[m2] += B3*cb[R3[m2]]; w[m0] += B1*cb[R1[m0]]
    r_ = w/v
    rhos[k] = (float(r_.max()), float(r_.max()/r_.min()-1))
    print(f"k={k}: rho={rhos[k][0]:.10f} cwgap={rhos[k][1]:.1e}", flush=True)

ks = np.array(sorted(rhos))
y = np.array([1 - rhos[k][0] for k in ks])
print("\n1-rho sequence:", " ".join(f"{x:.5f}" for x in y))
rat = y[1:]/y[:-1]
print("ratios:", " ".join(f"{r:.4f}" for r in rat))

# model 1: pure geometric on k>=8 (skip shallow transient)
sel = ks >= 8
lr = np.polyfit(ks[sel], np.log(y[sel]), 1)
r1 = np.exp(lr[0]); res1 = np.abs(np.log(y[sel]) - np.polyval(lr, ks[sel])).max()
print(f"\npure geometric fit (k>=8): rate={r1:.4f}, max log-resid={res1:.4f}")

# model 2: y_inf + C*r^k — grid over y_inf
best = None
for yinf in np.linspace(0, y[-1]*0.98, 200):
    z = y[sel] - yinf
    if (z <= 0).any(): continue
    lr2 = np.polyfit(ks[sel], np.log(z), 1)
    res = np.abs(np.log(z) - np.polyval(lr2, ks[sel])).max()
    if best is None or res < best[0]:
        best = (res, yinf, np.exp(lr2[0]))
print(f"leveling fit: y_inf={best[1]:.5f}, rate={best[2]:.4f}, max log-resid={best[0]:.4f}")
print(f"CEILING would need y_inf ~ 1-rho_inf(2) > 0; DENSITY: y_inf = 0")
print(f"note: ratios trend {'UP (leveling)' if rat[-1]>rat[len(rat)//2] else 'flat/down (geometric)'}")
