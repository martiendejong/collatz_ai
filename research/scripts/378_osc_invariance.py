# 378: Route A step 14 — oscillation invariance of the one-level tower map L.
# OSC_j(x) = RMS over 3-adic scale-j cells of std(log x within cell).
# If L maps an oscillation profile to a dominated profile (factor < 1 per
# scale, k-flat), the sub-cone {OSC_j <= C r^j} is invariant: the Doeblin
# ingredient. Inputs: (a) the real field cb3, (b) rough random fields.
import numpy as np
import os
from math import log2

ALPHA = log2(3.0)
CACHE = r"E:\projects\collatz\research\cache"

def osc_profile(x, jmax=5):
    lx = np.log(x)
    out = []
    for j in range(jmax):
        M = 3**(j+1)
        prof = lx.reshape(-1, M) if False else None
        cells = lx.reshape(len(lx)//M, M)
        out.append(float(np.sqrt((cells.std(axis=1)**2).mean())))
    return np.array(out)

for lam in [1.70, 2.00]:
    k = 12
    N = 3**(k-1); Nl = N//3; Nll = Nl//3
    s = np.arange(Nl, dtype=np.int64)
    sig4 = (4*s+2) % Nl
    sl = s//3
    R1p = (4*sl) % Nll; R3p = (2*sl+1) % Nll
    r = s % 3
    A = lam**-2; B1 = lam**(ALPHA-2); B3 = lam**(ALPHA-1)
    v = np.load(os.path.join(CACHE, f"v_lam{lam:.2f}_k{k}.npy")).astype(np.float64)
    rho = 1.0
    for _ in range(200):
        i = np.arange(N); T4 = (4*i+2) % N
        cbv = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        w = A*v[T4]
        ri = i % 3; si = i//3
        w[ri == 2] += B3*cbv[((2*si+1) % Nl)[ri == 2]]
        w[ri == 0] += B1*cbv[((4*si) % Nl)[ri == 0]]
        rho = float(w.max()); v = w/rho
    t = A/rho
    def L(h):
        f = np.zeros(Nl)
        f[r == 0] = (B1/rho)*h[R1p[r == 0]]
        f[r == 2] = (B3/rho)*h[R3p[r == 2]]
        g = np.zeros(Nl)
        for _ in range(int(np.log(1e-14)/np.log(t))+1):
            g = t*g[sig4] + f
        return g
    cbv = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
    cb3 = (cbv[:Nll] + cbv[Nll:2*Nll] + cbv[2*Nll:])/3.0
    rng = np.random.default_rng(378)
    rough = cb3*np.exp(0.5*rng.standard_normal(Nll))
    for name, h in [("echt cb3", cb3), ("ruw", rough)]:
        pin = osc_profile(h)
        pout = osc_profile(L(h))
        print(f"lam={lam} {name:>8}: osc-in  {np.array2string(pin, precision=3)}")
        print(f"{'':>14}  osc-uit {np.array2string(pout, precision=3)}  "
              f"ratio/schaal {np.array2string(pout/pin, precision=2)}", flush=True)
