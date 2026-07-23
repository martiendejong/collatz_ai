"""
171_renormalization_model.py
=============================
THE REDUCED RENORMALIZATION MODEL of the Krasikov hierarchy.

Exact structural input (Obs 318-321, all verified):
  - backbone x4 maps sikkel-depth a -> desert-depth a (4m+4 = 4(m+1));
    desert exits to D1 in one step (any depth); D1 re-enters the sikkel
    tower with geometric depth: P(v3(4m+1) = t) = (2/3) 3^(2-t), t >= 2.
  - D3 (sikkel depth a >= 2): c = L^-2 * c(desert) + L^(a-1-exponent) feed
    from depth a-1  [L := lambda, advanced coeff L^(alpha-1)]
  - D2 (desert):  c = L^-2 * c(D1)          [pure decay]
  - D1:           c = L^-2 * E[c(sikkel_t)] + L^(alpha-2) * E[c(feed_(t-1))]

Minimal state space: S_2..S_R (sikkel depths), D1, V (desert), with the
geometric re-entry distribution. Perron condition rho(lambda) = 1 gives
gamma_model(R) = log2 lambda*(R).

Questions:
  Q1: does gamma_model(R) land near the measured gamma(k) (0.84..0.902)?
  Q2: does gamma_model(R) -> 1 as R -> inf (the limit-theorem mechanism)?
  Q3: what is the R->inf closed-form value?
"""
import numpy as np
from math import log2

ALPHA = log2(3.0)

def rho(lmb, R, iters=4000):
    # state vector: [S_2..S_R, D1, V]
    n = (R - 1) + 2
    idxD1 = n - 2
    idxV = n - 1
    # geometric re-entry weights, capped at R (tail mass lumped on R)
    p = np.zeros(R + 1)
    for t in range(2, R):
        p[t] = (2/3) * 3.0**(2 - t)
    p[R] = 1.0 - p[2:R].sum()
    A2 = lmb**(-2.0)
    B3 = lmb**(ALPHA - 1.0)
    B1 = lmb**(ALPHA - 2.0)
    v = np.ones(n)
    g = 1.0
    for _ in range(iters):
        w = np.empty(n)
        # S_2 = A2*V + B3*D1 ; S_d = A2*V + B3*S_{d-1}
        w[0] = A2 * v[idxV] + B3 * v[idxD1]
        for d in range(3, R + 1):
            w[d - 2] = A2 * v[idxV] + B3 * v[d - 3]
        # V = A2 * D1
        w[idxV] = A2 * v[idxD1]
        # D1 = A2 * sum_t p_t S_t + B1 * (p_2 D1 + sum_{t>=3} p_t S_{t-1})
        ES = sum(p[t] * v[t - 2] for t in range(2, R + 1))
        EF = p[2] * v[idxD1] + sum(p[t] * v[t - 3] for t in range(3, R + 1))
        w[idxD1] = A2 * ES + B1 * EF
        g = np.abs(w).max()
        v = w / g
    return g

def gamma_model(R):
    lo, hi = 1.05, 1.999
    for _ in range(45):
        mid = 0.5 * (lo + hi)
        if rho(mid, R) >= 1.0:
            lo = mid
        else:
            hi = mid
    return log2(lo)

print("Reduced renormalization model: gamma_model(R)")
print(f"{'R':>4} {'gamma_model':>12}   vs measured gamma(k)")
meas = {11: 0.8417, 12: 0.8531, 13: 0.8630, 14: 0.8724, 15: 0.8812,
        16: 0.8893, 17: 0.8950, 18: 0.9020}
for R in [4, 6, 8, 10, 12, 14, 16, 18, 24, 32, 48, 64, 96, 128]:
    g = gamma_model(R)
    tag = f"  gamma({R}) = {meas[R]}" if R in meas else ""
    print(f"{R:>4} {g:>12.4f}{tag}")
print()
print("R -> inf behaviour: if gamma_model(R) -> 1, the reduced model already")
print("contains the limit-theorem mechanism; the remaining work is proving")
print("model fidelity (the 40% unexplained variance).")
