"""Prefix-increment transfer: the cascade recursion one level deeper.
m_p(X) = E[X|digits 0..p] - E[X|digits 0..p-1] (orthogonal increment fields).
Exact: rho*m_p(v) = A*m_p(Pv) + m_p(CBv).
Measure inc_p = Var(m_p) for v, Pv, CBv + cross-cov; and test the shift
hypothesis inc_p(Pv) ~ c * inc_{p+1}(v) (x4 pulls one digit down)."""
import numpy as np
from math import log2

ALPHA = log2(3.0)
CACHE = "E:/projects/collatz/research/cache"

for lam in [1.05, 1.70]:
    k = 13
    v = np.load(f"{CACHE}/v_lam{lam:.2f}_k{k}.npy")
    rho = float(open(f"{CACHE}/rho_lam{lam:.2f}_k{k}.txt").read())
    A = lam**-2; B1 = lam**(ALPHA-2); B3 = lam**(ALPHA-1)
    N = v.size; Nl = N//3
    i = np.arange(N, dtype=np.int64)
    T4 = (4*i+2) % N
    s_arr, r_arr = np.divmod(i, 3)
    R1 = (4*s_arr) % Nl; R3 = (2*s_arr+1) % Nl
    cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
    Pv = v[T4]
    CBv = np.zeros(N)
    CBv[r_arr==0] = B1*cb[R1[r_arr==0]]
    CBv[r_arr==2] = B3*cb[R3[r_arr==2]]

    P = k-1
    def inc_fields(X):
        """list of increment fields m_p (as reduced arrays of conditional means)"""
        Xc = X - X.mean()
        out = []
        prev = None
        for p in range(P):
            M = 3**(p+1)
            cm = Xc.reshape(N//M, M).mean(axis=0)      # E[X | low p+1 digits]
            if prev is None:
                m = cm.copy()
            else:
                m = cm - np.repeat(prev, 3)[:M]        # careful: prev has M/3 entries, tile pattern
            out.append((p, cm))
            prev = cm
        # convert to increments properly: m_p defined on M=3^{p+1} cells: cm_p - lift(cm_{p-1})
        incs = []
        prev = None
        for p, cm in out:
            if prev is None:
                m = cm
            else:
                # cell index c in [0,3^{p+1}): low digits; cm_{p-1} index = c mod 3^p
                M = cm.size
                m = cm - prev[np.arange(M) % (M//3)]
            incs.append(m)
            prev = cm
        return incs

    mv = inc_fields(v)
    mP = inc_fields(Pv)
    mC = inc_fields(CBv)

    print(f"\n=== lam={lam} k={k} ===")
    print("p | inc_p(v) | A^2*inc_p(Pv)/rho^2 | inc_p(CB)/rho^2 | 2A*cross/rho^2 | sum/inc_v")
    for p in range(P):
        iv_ = (mv[p]**2).mean()
        iP = (mP[p]**2).mean()
        iC = (mC[p]**2).mean()
        cr = (mP[p]*mC[p]).mean()
        tot = (A*A*iP + iC + 2*A*cr)/rho**2
        print(f"{p} | {iv_:.3e} | {A*A*iP/rho**2:.3e} | {iC/rho**2:.3e} | {2*A*cr/rho**2:+.3e} | {tot/iv_:.4f}")
    # shift hypothesis: inc_p(Pv) vs inc_{p+1}(v) and inc_p(v)
    print("shift test: inc_p(Pv)/inc_{p+1}(v) and /inc_p(v):")
    for p in range(P-1):
        iP = (mP[p]**2).mean()
        r1 = iP/((mv[p+1]**2).mean()+1e-300)
        r0 = iP/((mv[p]**2).mean()+1e-300)
        print(f"  p={p}: /inc_(p+1)={r1:.3f}  /inc_p={r0:.3f}")
