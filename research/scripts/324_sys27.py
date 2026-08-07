"""
Digit-2 fibered system (27 unknowns mu_r[u], u in Z/9) driven by 9 c-values:
  mu0[u] = t*mu2[4u] + (B1/rho)*c[4u]
  mu1[u] = t*mu0[4u+2]
  mu2[u] = t*mu1[4u+3] + (B3/rho)*c[2u+1]     (all maps mod 9)
Validate to 1e-15; catalog routing cycles; test: gap-richness by cb-orbit type
(6-cycle {0,1,3,7,6,4} vs 2-cycle {2,5} vs fixed {8}).
"""
import numpy as np
from math import log2
ALPHA = log2(3.0)
CACHE = "E:/projects/collatz/research/cache"

# cycle catalog
p4 = [(4*u) % 9 for u in range(9)]
p42 = [(4*u+2) % 9 for u in range(9)]
p43 = [(4*u+3) % 9 for u in range(9)]
p21 = [(2*u+1) % 9 for u in range(9)]
def cycles(p):
    seen = set(); out = []
    for s in range(9):
        if s in seen: continue
        c = [s]; seen.add(s); x = p[s]
        while x != s:
            c.append(x); seen.add(x); x = p[x]
        out.append(tuple(c))
    return out
print("cycle structures on Z/9:")
print("  4u   :", cycles(p4))
print("  4u+2 :", cycles(p42))
print("  4u+3 :", cycles(p43))
print("  2u+1 :", cycles(p21))

for lam in [1.05, 1.70]:
    k = 13
    v = np.load(f"{CACHE}/v_lam{lam:.2f}_k{k}.npy")
    rho = float(open(f"{CACHE}/rho_lam{lam:.2f}_k{k}.txt").read())
    A = lam**-2; B1 = lam**(ALPHA-2); B3 = lam**(ALPHA-1)
    t = A/rho
    N = v.size; Nl = N//3
    cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
    s = np.arange(Nl)
    su = s % 9
    c9 = np.array([cb[su == u].mean() for u in range(9)])
    mu_meas = np.zeros((3, 9))
    for r in range(3):
        vr = v[r::3]
        for u in range(9):
            mu_meas[r, u] = vr[su == u].mean()
    # solve the 27-system: mu0[u] = t^3*mu0[sigma3(u)] + f(u) with chained substitution
    # mu2[u] = t*mu1[4u+3] + (B3/rho)c[2u+1]; mu1[w] = t*mu0[4w+2]; so
    # mu0[u] = t*mu2[4u] + (B1/rho)c[4u]
    #        = t*(t*mu1[4*(4u)+3] + (B3/rho)c[2*(4u)+1]) + (B1/rho)c[4u]
    #        = t^2*mu0[4*(16u+3)+2] + t*(B3/rho)c[8u+1] + (B1/rho)c[4u]
    # sigma3(u) = (64u+14) mod 9 = (u + 5) mod 9   (64 mod 9 = 1, 14 mod 9 = 5)
    sig3 = [(u+5) % 9 for u in range(9)]
    f = np.array([t*(B3/rho)*c9[(8*u+1) % 9] + (B1/rho)*c9[(4*u) % 9] for u in range(9)])
    # mu0[u] = t^3 mu0[sig3(u)] + f(u): sig3 is u -> u+5 mod 9: a 9-cycle!
    # solve: mu0[u] = sum_{j=0}^{8} t^(3j) f(sig3^j(u)) / (1 - t^27)
    mu0 = np.zeros(9)
    for u in range(9):
        acc = 0.0; x = u
        for jj in range(9):
            acc += t**(3*jj) * f[x]
            x = sig3[x]
        mu0[u] = acc/(1 - t**27)
    mu1 = np.array([t*mu0[(4*u+2) % 9] for u in range(9)])
    mu2 = np.array([t*mu1[(4*u+3) % 9] + (B3/rho)*c9[(2*u+1) % 9] for u in range(9)])
    pred = np.stack([mu0, mu1, mu2])
    err = np.abs(pred - mu_meas).max()/mu_meas.mean()
    print(f"\nlam={lam}: 27-system max rel err = {err:.2e}")

    # gap-richness by cb-orbit type of u (digit-2 sub-class of the cb-index space)
    vbar = (v[:Nl]+v[Nl:2*Nl]+v[2*Nl:])/3.0
    G = vbar - cb
    rich = np.array([G[su == u].mean()/vbar[su == u].mean() for u in range(9)])
    six = [0,1,3,7,6,4]; two = [2,5]; fix = [8]
    print(f"  richness by orbit: 6-cycle {rich[six].mean():.5f} (sd {rich[six].std():.5f}) | "
          f"2-cycle {rich[two].mean():.5f} | fixed(8) {rich[fix].mean():.5f} | global {G.mean()/vbar.mean():.5f}")
