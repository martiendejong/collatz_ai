"""
Exact decomposition of the min-aggregation coefficient:
  m_q(cb) = m_q(vbar) - dG_q,  vbar = branch average (= m_q(v) cell-wise)
  regression dG = b*m_q(v) + eps (orthogonal)  =>  c_q = (1-b_q)^2 + Var(eps)/inc_q  EXACT.
Scale-invariance ansatz predicts eps ~ 0 and b_q -> const = relative top-gap coefficient.
"""
import numpy as np
from math import log2

ALPHA = log2(3.0)
CACHE = "E:/projects/collatz/research/cache"

for lam in [1.05, 1.30, 1.70]:
    k = 13
    v = np.load(f"{CACHE}/v_lam{lam:.2f}_k{k}.npy")
    N = v.size; Nl = N//3
    cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
    vbar3 = (v[:Nl] + v[Nl:2*Nl] + v[2*Nl:])/3.0    # branch average field on [0,Nl)
    Gfield = vbar3 - cb                               # pointwise top-triple gap >= 0

    P = k-2  # digits of the [0,Nl) index space
    def inc_fields(X, n):
        Xc = X - X.mean()
        out = []
        prev = None
        for p in range(P):
            M = 3**(p+1)
            cm = Xc.reshape(n//M, M).mean(axis=0)
            m = cm if prev is None else cm - prev[np.arange(M) % (M//3)]
            out.append(m)
            prev = cm
        return out

    m_vbar = inc_fields(vbar3, Nl)
    m_G = inc_fields(Gfield, Nl)
    m_cb = inc_fields(cb, Nl)

    print(f"\n=== lam={lam} k={k} ===")
    print("q | c_q | b_q | (1-b)^2 | Var(eps)/inc | check (1-b)^2+resid | corr(m,dG)")
    for q in range(P):
        mv = m_vbar[q]; mg = m_G[q]; mc = m_cb[q]
        inc_v = (mv**2).mean()
        c_q = (mc**2).mean()/inc_v
        b = (mv*mg).mean()/inc_v
        eps = mg - b*mv
        resid = (eps**2).mean()/inc_v
        chk = (1-b)**2 + resid
        corr = (mv*mg).mean()/np.sqrt((mv**2).mean()*(mg**2).mean()+1e-300)
        print(f"{q:2d} | {c_q:.4f} | {b:+.4f} | {(1-b)**2:.4f} | {resid:.4f} | {chk:.4f} | {corr:+.3f}")
    # global relative gap of top triples (the scale-proportionality constant candidate)
    print(f"global E[G]/E[vbar] = {Gfield.mean()/vbar3.mean():.5f}")
