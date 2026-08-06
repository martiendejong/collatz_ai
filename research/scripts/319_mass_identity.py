"""Exact mass identity: rho_lin - rho = ((B1+B3)/3) * E[G]/E[v], with
rho_lin = lam^-2 + (lam^(a-2)+lam^(a-1))/3 (uniform left eigenvector of the
mean-field operator; min<=mean gives rho <= rho_lin always).
Verify across cached (lam, k)."""
import numpy as np
from math import log2

ALPHA = log2(3.0)
CACHE = "E:/projects/collatz/research/cache"

print("lam | k | rho_lin | rho | deficit | ((B1+B3)/3)*E[G]/E[v] | ratio")
for lam in [1.05, 1.30, 1.70]:
    for k in [8, 12, 13, 14]:
        try:
            v = np.load(f"{CACHE}/v_lam{lam:.2f}_k{k}.npy")
            rho = float(open(f"{CACHE}/rho_lam{lam:.2f}_k{k}.txt").read())
        except Exception:
            continue
        A = lam**-2; B1 = lam**(ALPHA-2); B3 = lam**(ALPHA-1)
        rho_lin = A + (B1+B3)/3
        N = v.size; Nl = N//3
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        G = (v[:Nl]+v[Nl:2*Nl]+v[2*Nl:])/3 - cb
        rhs = (B1+B3)/3 * G.mean()/v.mean() * 3*Nl/N * 1.0  # E over j vs i: sum_j G / sum_i v = mean(G)*Nl/(mean(v)*N)
        # careful: deficit = (B1+B3) * sum_j G / sum_i v = (B1+B3)*mean(G)*Nl/(mean(v)*N) = (B1+B3)/3 * mean(G)/mean(v)
        rhs = (B1+B3)/3 * G.mean()/v.mean()
        deficit = rho_lin - rho
        print(f"{lam:.2f} | {k:2d} | {rho_lin:.6f} | {rho:.6f} | {deficit:.6f} | {rhs:.6f} | {deficit/rhs:.6f}")

# lam=2 series: deficit = 1 - rho(2,k) => E[G]/E[v] = (4/3)*(1-rho)
print("\nlam=2: 1-rho(2,k) = (3/4)*relgap  => relgap sequence:")
seq = [0.13899,0.12033,0.10498,0.08967,0.08061,0.07327,0.06651,0.06016,0.05495,0.05002,0.04550,0.04148]
for k, y in zip(range(5,17), seq):
    print(f"  k={k}: relgap = {y*4/3:.5f}")
