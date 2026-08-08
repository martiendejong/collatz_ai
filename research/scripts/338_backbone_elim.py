"""
Backbone elimination theorem (candidate): from rho*v = A*v(T4) + feed(cb),
iterate the backbone: v(i) = (1/rho) * sum_{n>=0} t^n * feed(cb)(T4^n(i)).
Since T4 advances class 0->2->1 with exact period 3, the feed-coefficient pattern
along the single T4-cycle is 3-periodic. Verify the Neumann reconstruction to
machine precision, then the cb-only fixed-point equation w = min3(K w).
"""
import numpy as np
from math import log2
CACHE = "E:/projects/collatz/research/cache"
ALPHA = log2(3.0)

for lam, k in [(1.05, 10), (1.70, 10)]:
    # need a k=10 vector: compute quickly
    A = lam**-2; B1 = lam**(ALPHA-2); B3 = lam**(ALPHA-1)
    N = 3**(k-1); Nl = N//3
    i = np.arange(N, dtype=np.int64)
    T4 = (4*i+2) % N
    s_arr, r_arr = np.divmod(i, 3)
    m0 = r_arr==0; m2 = r_arr==2
    R1 = (4*s_arr) % Nl; R3 = (2*s_arr+1) % Nl
    v = np.ones(N); rho = 1.0
    for _ in range(1200):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        w = A*v[T4]; w[m2] += B3*cb[R3[m2]]; w[m0] += B1*cb[R1[m0]]
        rho = float(w.max()); w /= rho; v = w
    t = A/rho
    cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])

    # feed field f(i): B1*cb[R1[s]] if i=3s; B3*cb[R3[s]] if i=3s+2; 0 if class1
    f = np.zeros(N)
    f[m0] = B1*cb[R1[m0]]
    f[m2] = B3*cb[R3[m2]]

    # Neumann reconstruction: v_rec(i) = (1/rho) sum_n t^n f(T4^n(i)), truncate when t^n < 1e-18
    nmax = int(np.ceil(np.log(1e-18)/np.log(t)))
    v_rec = np.zeros(N)
    idx = i.copy()
    coef = 1.0/rho
    for n in range(nmax):
        v_rec += coef * f[idx]
        idx = T4[idx]
        coef *= t
    err = np.abs(v_rec - v).max()/v.mean()
    # cycle check
    # verify T4 single-cycle claim quickly on this N
    seen = np.zeros(N, dtype=bool)
    x = 0; L = 0
    while not seen[x]:
        seen[x] = True; x = int(T4[x]); L += 1
    print(f"lam={lam} k={k}: Neumann-reconstructie max rel err = {err:.2e} "
          f"(nmax={nmax}) | T4-cykel vanaf 0 heeft lengte {L} van N={N} "
          f"{'(SINGLE CYCLE)' if L==N else '(meerdere cykels!)'}")
