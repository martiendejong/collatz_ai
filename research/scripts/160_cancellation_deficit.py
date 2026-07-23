"""
160_cancellation_deficit.py
============================
THE CANCELLATION DEFICIT -- the mastermind panel's prescribed measurement.

Exact requirement for a Fourier-side proof of N(0) = 0 at signature (A,B):
    |sum_{t!=0} S(t)| < D - S(0)                                     (+)
(then N(0) = (1/D)|S(0) + sum_{t!=0}S(t)| < 1, and N(0) integer => 0).

In reality the cancellation is PERFECT: sum_{t!=0} S(t) = -S(0) exactly
(since N(0)=0). The question is which BOUND could ever prove (+):

  L1 route:      sum |S(t)| < D - S(0)?    -> measure L1 = sum_{t!=0}|S(t)|
  L2 benchmark:  Parseval: sum_t |S(t)|^2 = D * sum_r N(r)^2  (exact check)
                 Cauchy-Schwarz upper bd for L1: sqrt((D-1)(L2 - S(0)^2))
  concentration: how many t carry 50% / 90% of the L1 mass? are they
                 structurally identifiable ("major arcs": small doubling-orbit
                 representatives / small t)?

Verdict metric: deficit = L1 / (D - S(0)). If deficit >> 1, NO absolute-value
bound (circle method with |.|) can ever prove (+): a proof must capture
sign-level (phase) cancellation -- and by how much (the deficit factor).
"""
import numpy as np
from math import comb, log2, gcd

THETA = log2(3)

def N_vector(A, B, M):
    w = [pow(3, A-1-i, M) for i in range(A)]
    dp = [np.zeros(M, dtype=np.int32) for _ in range(A+1)]
    dp[1][w[0] % M] = 1
    for d in range(1, B):
        p2 = pow(2, d, M)
        for i in range(min(d, A-1), 0, -1):
            dp[i+1] += np.roll(dp[i], (w[i]*p2) % M)
    return dp[A]

for A, B in [(10, 16), (12, 20), (13, 21)]:
    D = (1 << B) - 3**A
    N = N_vector(A, B, D)
    S0 = comb(B-1, A-1)
    S = np.fft.ifft(N.astype(np.float64)) * D
    absS = np.abs(S)
    L1 = float(absS[1:].sum())
    L2 = float((absS**2).sum())
    parseval = D * float((N.astype(np.float64)**2).sum())
    need = D - S0
    cs = np.sqrt((D-1) * (L2 - S0**2))
    print("="*74)
    print(f"(A,B)=({A},{B})  D={D:,}  S(0)={S0:,}  threshold D-S(0)={need:,}")
    print("="*74)
    print(f"Parseval check: sum|S|^2 = {L2:.6g} vs D*sum N^2 = {parseval:.6g}  "
          f"(match: {abs(L2-parseval)/parseval < 1e-9})")
    print(f"actual sum_(t!=0) S(t) = {float(S[1:].sum().real):.2f}  "
          f"(exact requirement: -S(0) = {-S0})")
    print(f"L1 mass  sum|S(t)|      = {L1:,.0f}")
    print(f"CANCELLATION DEFICIT: L1/(D-S(0)) = {L1/need:,.1f}x")
    print(f"Cauchy-Schwarz would give {cs:,.0f}  ({cs/need:,.1f}x threshold)")
    # concentration of L1 mass
    order_idx = np.argsort(absS[1:])[::-1] + 1
    csum = np.cumsum(absS[order_idx])
    n50 = int(np.searchsorted(csum, 0.5*L1)) + 1
    n90 = int(np.searchsorted(csum, 0.9*L1)) + 1
    print(f"L1 concentration: top {n50:,} of {D-1:,} frequencies carry 50%; "
          f"top {n90:,} carry 90%")
    # are top frequencies structurally simple? (small representative in <2>-orbit)
    def min_in_orbit(t):
        best = t; x = t
        for _ in range(200):
            x = 2*x % D
            m = min(x, D-x)
            if m < best: best = m
            if x == t: break
        return best
    top = order_idx[:10]
    print("top-10 |S| frequencies: t, |S(t)|/S(0), gcd(t,D), min|.| in 2-orbit:")
    for t in top:
        print(f"   t={int(t):>7} |S|/S0={absS[t]/S0:.4f} gcd={gcd(int(t),D):>6} "
              f"orbit-min={min_in_orbit(int(t)):>7}")
    print()
