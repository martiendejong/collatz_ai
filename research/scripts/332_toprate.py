"""Precision instrument: top-layer across-k rate = inc_{k-2}(k+1)/inc_{k-2->}(k)...
define TR(k) = inc_last(k+1)/inc_last(k). Measure at all lambdas with k=12,13,14
(compute missing vectors), then re-mine closed forms on the clean c(lam) values."""
import numpy as np
from math import log2
import os
CACHE = "E:/projects/collatz/research/cache"
ALPHA = log2(3.0)

def get_v(lam, k):
    fn = f"{CACHE}/v_lam{lam:.2f}_k{k}.npy"
    if os.path.exists(fn):
        return np.load(fn)
    A = lam**-2; B1 = lam**(ALPHA-2); B3 = lam**(ALPHA-1)
    N = 3**(k-1); Nl = N//3
    i = np.arange(N, dtype=np.int64)
    T4 = (4*i+2) % N
    s_arr, r_arr = np.divmod(i, 3)
    m0 = r_arr==0; m2 = r_arr==2
    R1 = (4*s_arr) % Nl; R3 = (2*s_arr+1) % Nl
    v = np.ones(N)
    for _ in range(500 if k <= 12 else 350):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        w = A*v[T4]; w[m2] += B3*cb[R3[m2]]; w[m0] += B1*cb[R1[m0]]
        v = w/w.max()
    np.save(fn, v)
    return v

def inc_last(lam, k):
    v = get_v(lam, k)
    N = v.size
    F = np.log2(v); F -= F.mean()
    prev = None; last = None
    for p in range(k-1):
        M = 3**(p+1)
        cm = F.reshape(N//M, M).mean(axis=0)
        m = cm if prev is None else cm - prev[np.arange(M) % (M//3)]
        last = float((m**2).mean()); prev = cm
    return last

LAMS = [1.05, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 1.95, 2.00]
rows = []
for lam in LAMS:
    i12 = inc_last(lam, 12); i13 = inc_last(lam, 13); i14 = inc_last(lam, 14)
    tr1 = i13/i12; tr2 = i14/i13
    rows.append((lam, tr1, tr2))
    print(f"lam={lam:.2f}: TR(12->13)={tr1:.4f}  TR(13->14)={tr2:.4f}  drift={tr2-tr1:+.4f}", flush=True)

# closed-form mining on TR(13->14) as clean c-estimate
print("\nmining c(lam) = TR against candidates (train even idx, test odd):")
import numpy as np
data = [(lam, tr2) for lam, _, tr2 in rows]
# need t, rho per lam at k=14
def eig(lam, k=14):
    fr = f"{CACHE}/rho_lam{lam:.2f}_k{k}.txt"
    if os.path.exists(fr):
        return float(open(fr).read())
    # quick CW from cached v
    v = get_v(lam, k)
    A = lam**-2; B1 = lam**(ALPHA-2); B3 = lam**(ALPHA-1)
    N = v.size; Nl = N//3
    i = np.arange(N, dtype=np.int64)
    T4 = (4*i+2) % N
    s_arr, r_arr = np.divmod(i, 3)
    m0 = r_arr==0; m2 = r_arr==2
    R1 = (4*s_arr) % Nl; R3 = (2*s_arr+1) % Nl
    cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
    w = A*v[T4]; w[m2] += B3*cb[R3[m2]]; w[m0] += B1*cb[R1[m0]]
    rho = float((w/v).max())
    open(fr, 'w').write(repr(rho))
    return rho

full = []
for lam, tr in data:
    rho = eig(lam)
    t = lam**-2/rho
    full.append(dict(lam=lam, c=tr, t=t, rho=rho))
tr_set = full[0::2]; te_set = full[1::2]
cands = [
  "t + (1-t)/3", "(1+2*t)/3", "(1+t)/2", "1-(1-t)**2/2",
  "(2+t)/3 - t*t/3", "2/3 + t/3", "1/(1+ (1-t)**2)", "1 - (1-t)*(1-t/2)",
  "(1+t+t*t)/ (1+t+1)", "3/(4-t)", "1/(2-t)", "(1+t)/(2-t*t)", "t**0.25",
]
for expr in cands:
    try:
        etr = max(abs(eval(expr, {}, r)-r['c'])/r['c'] for r in tr_set)
        ete = max(abs(eval(expr, {}, r)-r['c'])/r['c'] for r in te_set)
        tag = "  <<< CANDIDATE" if (etr < 0.01 and ete < 0.01) else ""
        if etr < 0.03:
            print(f"  c = {expr}: {etr:.4f} | {ete:.4f}{tag}")
    except Exception:
        pass
