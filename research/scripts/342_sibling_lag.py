"""
Sibling-lag theorem (candidate):
(1) Translation by Nl commutes with T4 (algebra: 4*Nl = Nl mod N) => the cycle
    lag between i and i+Nl is a GLOBAL constant L0 with 3*L0 = 0 mod N, so
    L0 in {N/3, 2N/3}, giving v3(L0) = k-2: MAXIMAL valuation.
(2) Hence sibling correlation = ACF(L0) = phi(k-2) = top of the valuation ladder
    -> explains rho_intra -> 1 (Obs 481) exactly.
Verify numerically at k=10, 12.
"""
import numpy as np
from math import log2
CACHE = "E:/projects/collatz/research/cache"

for lam, k in [(1.05, 12), (1.70, 12), (2.00, 12)]:
    v = np.load(f"{CACHE}/v_lam{lam:.2f}_k{k}.npy")
    N = v.size; Nl = N//3
    i = np.arange(N, dtype=np.int64)
    T4 = (4*i+2) % N
    # commutation check (algebra says exact)
    comm = np.abs(T4[(i+Nl) % N] - (T4 + Nl) % N).max()
    # find L0: walk the cycle from 0 until reaching Nl
    x = 0; L0 = 0
    target = Nl
    while True:
        x = int(T4[x]); L0 += 1
        if x == (0 + target) % N:
            break
        if L0 > N:
            L0 = -1; break
    # sibling correlation vs ACF at that lag
    F = np.log2(v); F -= F.mean()
    var = float((F**2).mean())
    sib = float((F*np.roll(F, -Nl)).mean())/var  # corr(F(i), F(i+Nl))
    # v3 of L0
    def v3(n):
        c = 0
        while n % 3 == 0:
            n //= 3; c += 1
        return c
    print(f"lam={lam} k={k}: commutatie-check={comm} | L0={L0} (N/3={Nl}) "
          f"v3(L0)={v3(L0) if L0>0 else '?'} (k-2={k-2}) | sibling-corr={sib:+.4f}")
