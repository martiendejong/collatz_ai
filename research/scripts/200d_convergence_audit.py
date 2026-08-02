"""
200d_convergence_audit.py
=========================
Self-audit of the g-column before any analytic work: is the +0.003
increment plateau (Obs 403) real, or an iteration-convergence artifact?

All g-points were computed with 250-300 power iterations. Growth
converges to 8 digits well before that, but the FINE-SCALE statistic
Var_end (exactly what V_k measures) can equilibrate more slowly -- and
the plateau is only 0.4% per depth step. If Var_end(k, 1.70) drifts by
more than ~0.1% between 300 and 1200 iterations, the plateau conclusion
is compromised and every g-point needs re-measurement at tight
convergence.

Measured here: Var_end(k, 1.70) at n_iter in {150, 300, 600, 1200} for
k = 15, 16, 17, plus the resulting d_15/d_16 per iteration budget.
"""
import numpy as np
from math import log2

ALPHA = log2(3.0)
LAM = 1.70
A, B1, B3 = LAM ** -2.0, LAM ** (ALPHA - 2.0), LAM ** (ALPHA - 1.0)


def maps(k):
    N = 3 ** (k - 1)
    i = np.arange(N, dtype=np.int64)
    T4 = (4 * i + 2) % N
    s, r = np.divmod(i, 3)
    Nl = N // 3
    R1m = ((4 * s) % Nl)[r == 0]
    R3m = ((2 * s + 1) % Nl)[r == 2]
    return N, Nl, T4, (r == 0), R1m, (r == 2), R3m


def var_end_at(k, iters_list):
    N, Nl, T4, m1, R1m, m3, R3m = maps(k)
    v = np.ones(N, dtype=np.float64)
    out = {}
    done = 0
    for target in iters_list:
        for _ in range(target - done):
            cb = np.minimum(np.minimum(v[:Nl], v[Nl:2 * Nl]), v[2 * Nl:])
            w2 = A * v[T4]
            w2[m1] += B1 * cb[R1m]
            w2[m3] += B3 * cb[R3m]
            v = w2 / w2.max()
        done = target
        T = np.stack([v[:Nl], v[Nl:2 * Nl], v[2 * Nl:]])
        X = np.log2(T) - np.log2(T.mean(axis=0))[None, :]
        out[target] = float(np.var(X))
    return out


ITS = (150, 300, 600, 1200)
res = {}
for k in (15, 16, 17):
    res[k] = var_end_at(k, ITS)
    print(f"  k={k}: " + "  ".join(f"it{t}={res[k][t]:.8f}" for t in ITS),
          flush=True)
    base = res[k][1200]
    print(f"        rel drift vs it1200: " +
          "  ".join(f"{(res[k][t]/base-1)*100:+.4f}%" for t in ITS),
          flush=True)
for t in ITS:
    d15 = res[16][t] / res[15][t]
    d16 = res[17][t] / res[16][t]
    print(f"  it{t}: d_15={d15:.5f}  d_16={d16:.5f}", flush=True)
print("done", flush=True)
