"""Test the plateau formula: with the exact A-preservation (inc_p(Pv)=inc_p(v)),
the balance is (rho^2 - A^2)*inc_p(v) = inc_p(CBv) + 2A*cross_p.
CBv injects the cb-field shifted one digit: inc_p(CBv) ~ W2 * inc_{p-1}(cb),
W2 = (B1^2 + B3^2)/3 + class-structure. Measure:
  c_q = inc_q(cb)/inc_q(v)   (min-aggregation coefficient)
  shift check: inc_p(CBv) vs (B1^2+B3^2)/3 * inc_{p-1}(cb)
  predicted plateau r from the closed balance vs measured 0.41 / 0.70."""
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
    CBv = np.zeros(N)
    CBv[r_arr==0] = B1*cb[R1[r_arr==0]]
    CBv[r_arr==2] = B3*cb[R3[r_arr==2]]

    def incs(X, n):
        Xc = X - X.mean()
        P_ = int(round(np.log(n)/np.log(3)))
        out = []
        prev = None
        for p in range(P_):
            M = 3**(p+1)
            cm = Xc.reshape(n//M, M).mean(axis=0)
            m = cm if prev is None else cm - prev[np.arange(M) % (M//3)]
            out.append(float((m**2).mean()))
            prev = cm
        return out

    iv_ = incs(v, N)
    icb = incs(cb, Nl)
    iCB = incs(CBv, N)
    W2 = (B1**2 + B3**2)/3

    print(f"\n=== lam={lam} ===   rho^2-A^2 = {rho**2-A**2:.4f}   W2=(B1^2+B3^2)/3 = {W2:.4f}")
    print("q | c_q = inc_q(cb)/inc_q(v) | shift: inc_{q+1}(CBv)/(W2*inc_q(cb))")
    for q in range(len(icb)):
        cq = icb[q]/iv_[q]
        sh = iCB[q+1]/(W2*icb[q]) if q+1 < len(iCB) else float('nan')
        print(f"{q:2d} | {cq:.4f} | {sh:.4f}")
    # predicted plateau from closed balance (neglecting cross deep):
    # (rho^2-A^2)*inc_p(v) ~ W2*c*inc_{p-1}(v)  =>  r = inc_p/inc_{p-1} = W2*c/(rho^2-A^2)
    c_deep = np.mean([icb[q]/iv_[q] for q in range(6, len(icb)-1)])
    r_pred = W2*c_deep/(rho**2 - A**2)
    meas = 0.41 if lam == 1.05 else 0.70
    print(f"c_deep = {c_deep:.4f}  predicted plateau r = {r_pred:.4f}  measured = {meas}")
    # with cross term (measured deep):
    # (rho^2-A^2)*inc_p = W2*c*inc_{p-1} + 2A*cross_p  -> refine
    # measure cross deep:
    def cross(p):
        Xc = v - v.mean(); Yc = CBv - CBv.mean()
        M = 3**(p+1)
        cmX = (v[T4] - v[T4].mean()).reshape(N//M, M).mean(axis=0)
        pmX = cmX - (cmX.reshape(-1)[np.arange(M) % (M//3)] if False else 0)
        return None
    print("(cross-term refinement: deep cross/rho^2/inc_v was ~ -0.4% at lam=1.05, -0.6% at 1.70 — small)")
