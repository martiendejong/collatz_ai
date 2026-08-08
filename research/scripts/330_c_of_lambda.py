"""
Systematic: measure the cascade plateau r(lambda) and implied c(lambda) across the
full lambda range at k=13, then relation-mine closed-form candidates with
train/test discipline (train: even-indexed lambdas, test: odd-indexed).
Edge rule: plateau = median of inc-ratios at digits p=4..7 (bulk only).
"""
import numpy as np
from math import log2
import os
CACHE = "E:/projects/collatz/research/cache"
ALPHA = log2(3.0)

def get_vr(lam, k=13):
    fn = f"{CACHE}/v_lam{lam:.2f}_k{k}.npy"
    fr = f"{CACHE}/rho_lam{lam:.2f}_k{k}.txt"
    A = lam**-2; B1 = lam**(ALPHA-2); B3 = lam**(ALPHA-1)
    if os.path.exists(fn) and os.path.exists(fr):
        return np.load(fn), float(open(fr).read())
    N = 3**(k-1); Nl = N//3
    i = np.arange(N, dtype=np.int64)
    T4 = (4*i+2) % N
    s_arr, r_arr = np.divmod(i, 3)
    m0 = r_arr==0; m2 = r_arr==2
    R1 = (4*s_arr) % Nl; R3 = (2*s_arr+1) % Nl
    v = np.ones(N); rho = 1.0
    for _ in range(450):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        w = A*v[T4]; w[m2] += B3*cb[R3[m2]]; w[m0] += B1*cb[R1[m0]]
        rho = float(w.max()); w /= rho; v = w
    np.save(fn, v); open(fr, 'w').write(repr(rho))
    return v, rho

def plateau(v, k=13):
    N = v.size
    F = np.log2(v); F -= F.mean()
    incs = []; prev = None
    for p in range(k-1):
        M = 3**(p+1)
        cm = F.reshape(N//M, M).mean(axis=0)
        m = cm if prev is None else cm - prev[np.arange(M) % (M//3)]
        incs.append(float((m**2).mean())); prev = cm
    rats = [incs[p+1]/incs[p] for p in range(len(incs)-1)]
    return float(np.median(rats[4:8]))   # bulk digits, edge-excluded

rows = []
LAMS = [1.05, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 1.95, 2.00]
for lam in LAMS:
    v, rho = get_vr(lam)
    A = lam**-2; B1 = lam**(ALPHA-2); B3 = lam**(ALPHA-1)
    t = A/rho
    W2 = (B1**2 + B3**2)/3
    r = plateau(v)
    c = r*(rho**2 - A**2)/W2
    N = v.size; Nl = N//3
    cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
    gam = 1 - cb.mean()/v.mean()*3*Nl/N*(v.mean()*N)/(v.mean()*N)  # = 1 - c̄/mean over [0,Nl)? fix:
    vbar = (v[:Nl]+v[Nl:2*Nl]+v[2*Nl:])/3.0
    gam = float((vbar-cb).mean()/vbar.mean())
    R = (t*t+lam)/(1+t*lam)
    rows.append(dict(lam=lam, t=t, rho=rho, r=r, c=c, gam=gam, R=R, W2=W2))
    print(f"lam={lam:.2f}: t={t:.4f} rho={rho:.5f} r={r:.4f} c={c:.4f} gam={gam:.5f}", flush=True)

# relation mining with train/test
import itertools
tr = rows[0::2]; te = rows[1::2]
def ev(expr, row):
    lam, t, R, gam, rho = row['lam'], row['t'], row['R'], row['gam'], row['rho']
    return eval(expr)
cands = [
  "1 - gam", "(1-gam)**2", "1 - 2*gam", "1/(1+2*gam)",
  "1 - gam*(1+lam)", "1 - lam*gam", "(1-gam)**lam",
  "1 - gam/(1-t)", "1 - gam*(1+t)/(1-t)",
  "1 - 3*gam", "1 - gam*(2+t)", "1 - gam*lam/(lam-1) if lam>1.01 else 0",
  "1 - gam*(1+1/t)", "1-2*lam*gam/(1+t)",
]
print("\nrelation mining for c(lam) [train err | test err], only <2% shown:")
for expr in cands:
    try:
        etr = max(abs(ev(expr, rw)-rw['c'])/rw['c'] for rw in tr)
        ete = max(abs(ev(expr, rw)-rw['c'])/rw['c'] for rw in te)
        if etr < 0.02 and ete < 0.02:
            print(f"  c = {expr}:  {etr:.4f} | {ete:.4f}   <<< CANDIDATE")
        elif etr < 0.05:
            print(f"  c = {expr}:  {etr:.4f} | {ete:.4f}")
    except Exception:
        pass
