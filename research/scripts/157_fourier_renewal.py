"""
157_fourier_renewal.py
=======================
THE FOURIER SIDE OF THE AFFINE RENEWAL IDENTITY.

Characters: S(t) = sum_r N(r) e^(2 pi i t r / D)  (so N(0) = (1/D) sum_t S(t)).
The renewal bijection C(sigma E) = 2C(E) - r*, r* = 3^(A-1), gives EXACTLY:

    e^(-2 pi i t r*/D) [ S(2t) - S_last(2t) ]  =  S(t) - S_first(t)      (**)

i.e. S evolves along the DOUBLING orbit t -> 2t mod D, driven by boundary
character sums. Since N(0) = 0 at computable signatures, we have the exact
cancellation sum_{t != 0} S(t) = -S(0) = -#patterns. Questions:

 1. Verify (**) numerically (float precision) at random t.
 2. Decompose {1..D-1} into doubling orbits; compute each orbit's sum of S(t).
    Is the mass of the cancellation concentrated on few orbits (structure to
    exploit) or spread across all of them (another honest wall)?
"""
import numpy as np
from math import gcd, comb, log2

THETA = log2(3)

def vectors(A, B, M):
    w = [pow(3, A-1-i, M) for i in range(A)]
    dp = [np.zeros(M, dtype=np.int64) for _ in range(A+1)]
    dp[1][w[0] % M] = 1
    dpF = [np.zeros(M, dtype=np.int64) for _ in range(A+1)]
    for d in range(1, B):
        p2 = pow(2, d, M)
        for i in range(min(d, A-1), 0, -1):
            c = (w[i]*p2) % M
            dp[i+1] += np.roll(dp[i], c)
            if i >= 2:
                dpF[i+1] += np.roll(dpF[i], c)
            if i == 1 and d == 1:
                dpF[2] += np.roll(dp[1], c)
    dp2 = [np.zeros(M, dtype=np.int64) for _ in range(A+1)]
    dp2[1][w[0] % M] = 1
    for d in range(1, B-1):
        p2 = pow(2, d, M)
        for i in range(min(d, A-1), 0, -1):
            dp2[i+1] += np.roll(dp2[i], (w[i]*p2) % M)
    NL = np.roll(dp2[A-1], (w[A-1]*pow(2, B-1, M)) % M)
    return dp[A], dpF[A], NL

for A, B in [(10, 16), (12, 20)]:
    D = (1 << B) - 3**A
    N, NF, NL = vectors(A, B, D)
    S  = np.fft.ifft(N.astype(np.float64))  * D   # S(t) = sum N(r) e^{+2pi i tr/D}
    SF = np.fft.ifft(NF.astype(np.float64)) * D
    SL = np.fft.ifft(NL.astype(np.float64)) * D
    rstar = pow(3, A-1, D)

    print("="*74)
    print(f"(A,B)=({A},{B})  D={D}  S(0)={S[0].real:.0f}  N(0)={int(N[0])}")
    print("="*74)

    # 1. verify (**) at random t
    rng = np.random.default_rng(0)
    ts = rng.integers(1, D, 12)
    worst = 0.0
    for t in ts:
        t = int(t); t2 = (2*t) % D
        lhs = np.exp(-2j*np.pi*t*rstar/D) * (S[t2] - SL[t2])
        rhs = S[t] - SF[t]
        worst = max(worst, abs(lhs - rhs))
    print(f"character-form renewal identity (**): max |lhs-rhs| over 12 random t "
          f"= {worst:.2e}  (float-exact)")

    # 2. doubling-orbit decomposition and orbit sums
    visited = np.zeros(D, dtype=bool)
    visited[0] = True
    orbit_sums = []
    for t0 in range(1, D):
        if visited[t0]:
            continue
        t = t0; osum = 0.0+0.0j; length = 0
        while not visited[t]:
            visited[t] = True
            osum += S[t]; length += 1
            t = (2*t) % D
        orbit_sums.append((abs(osum), osum.real, length, t0))
    orbit_sums.sort(reverse=True)
    total = sum(o[1] for o in orbit_sums)
    print(f"number of doubling orbits: {len(orbit_sums)}; "
          f"sum of all orbit sums = {total:.2f} (must be -S(0) = {-S[0].real:.0f})")
    print(f"top orbit |sums| (abs, real, length, representative):")
    for a, re, ln, t0 in orbit_sums[:8]:
        print(f"   |sum|={a:12.1f}  re={re:12.1f}  len={ln:>6}  t0={t0}"
              f"   gcd(t0,D)={gcd(t0, D)}")
    # concentration measure
    top1 = orbit_sums[0][0]
    frac = top1 / sum(a for a, _, _, _ in orbit_sums)
    print(f"concentration: largest orbit carries {100*frac:.1f}% of total |mass|")
    print()
