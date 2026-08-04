"""
207c_d3_rate.py
===============
Obs 412 follow-up: instead of raw D3 count or log-weight (both conflate
D3 signal with walk length), use the D3 RATE = #D3 / alive_steps.

This normalises out walk-length variation and should give a cleaner
common-cause proxy.  The walk terminates when D2 is encountered (alive
flag goes false); among the live steps, the D3 fraction is the quantity
that drives both level and spread monotonically.

Prediction: partial_corr(logL, logS | D3_rate) <= partial with raw D3
count (0.66), i.e. the rate explains at least as much of the coupling.
"""
import numpy as np
from math import log2

ALPHA = log2(3.0)
LAM   = 1.70
A     = LAM ** -2.0
B1    = LAM ** (ALPHA - 2.0)
B3    = LAM ** (ALPHA - 1.0)
G     = 8


def perron_and_d3rate(k, n_iter=300):
    N  = 3 ** (k - 1)
    i  = np.arange(N, dtype=np.int64)
    T4 = (4 * i + 2) % N
    s, r = np.divmod(i, 3)
    Nl   = N // 3
    R1   = (4 * s) % Nl
    R3   = (2 * s + 1) % Nl
    m1, m3 = (r == 0), (r == 2)

    v = np.ones(N, dtype=np.float64)
    for _ in range(n_iter):
        cb   = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        w2   = A * v[T4]
        w2[m1] += B1 * cb[R1[m1]]
        w2[m3] += B3 * cb[R3[m3]]
        v    = w2 / w2.max()
    v = v / v.mean()

    pos       = np.arange(Nl, dtype=np.int64)
    d3count   = np.zeros(Nl, dtype=np.float64)
    alive_cnt = np.zeros(Nl, dtype=np.float64)
    alive     = np.ones(Nl, dtype=bool)
    mod       = Nl
    for _ in range(G):
        t   = pos % 3
        alive_cnt[alive] += 1.0
        d3count[(t == 2) & alive] += 1.0
        alive = alive & (t != 1)
        sp  = pos // 3
        mod //= 3
        nxt = np.where(t == 0, (4 * sp) % mod, (2 * sp + 1) % mod)
        pos = np.where(alive, nxt, pos)

    # rate: #D3 / alive_steps; 0 if no live steps
    d3rate = np.where(alive_cnt > 0, d3count / alive_cnt, 0.0)
    return v, d3rate


print(f"D3-rate common-cause test  lam={LAM}  G={G}", flush=True)
print("  sys  corr(L,S)  partial|rate   partial_raw", flush=True)

PARTIAL_RAW = {15: 0.65995, 16: 0.65851, 17: 0.65776}

for kp1 in (15, 16, 17):
    v, d3rate = perron_and_d3rate(kp1)
    Nl  = v.size // 3
    mn  = (v[:Nl] + v[Nl:2*Nl] + v[2*Nl:]) / 3.0
    Nl3 = Nl // 3
    M   = np.stack([mn[:Nl3], mn[Nl3:2*Nl3], mn[2*Nl3:]])
    level  = M.mean(axis=0)
    spread = M.std(axis=0)
    W      = np.stack([d3rate[:Nl3], d3rate[Nl3:2*Nl3], d3rate[2*Nl3:]]).mean(axis=0)

    ok  = spread > 0
    x   = np.log(level[ok]);  x -= x.mean()
    y   = np.log(spread[ok]); y -= y.mean()
    z   = W[ok].astype(np.float64); z -= z.mean()

    def c(a, b):
        return float(np.mean(a*b) / np.sqrt(np.mean(a*a)*np.mean(b*b)))

    cxy, cxz, cyz = c(x, y), c(x, z), c(y, z)
    partial = (cxy - cxz*cyz) / np.sqrt((1 - cxz**2)*(1 - cyz**2))

    print(f"  {kp1:2d}  {cxy:.5f}    {partial:.5f}        "
          f"{PARTIAL_RAW[kp1]:.5f}", flush=True)

print("done", flush=True)
