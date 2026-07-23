"""
158_fourier_structure.py
=========================
FOLLOW-UP on the identical orbit sums found in script 157.

 Q1. How many doubling-orbits share the identical sum? Which symmetry group
     explains the multiplicity? Candidates: t -> -t (conjugation; -1 in <2>?)
     and t -> 3t (note 3^A == 2^B mod D, so the image of 3 in the quotient
     Q = (Z/D)*/<2> has order dividing A).
 Q2. Is N itself invariant under r -> 3r mod D (which would force S(3t)=S(t))?
 Q3. Divisor-level counts: Z_d = #{patterns : d | C} for every divisor d | D,
     vs the random model S/d. The CRT view: all Fourier mass at gcd(t,D)=D/d
     is determined by N mod d. Deviations Z_d - S/d measure arithmetic bias
     prime-by-prime -- the entire cycle problem is Z_D = 0 vs S/D ~ 2^-0.05B.
"""
import numpy as np
from math import gcd, comb, log2

THETA = log2(3)

def vectors_N(A, B, M):
    w = [pow(3, A-1-i, M) for i in range(A)]
    dp = [np.zeros(M, dtype=np.int64) for _ in range(A+1)]
    dp[1][w[0] % M] = 1
    for d in range(1, B):
        p2 = pow(2, d, M)
        for i in range(min(d, A-1), 0, -1):
            dp[i+1] += np.roll(dp[i], (w[i]*p2) % M)
    return dp[A]

def divisors(n, fac):
    ds = [1]
    for p, e in fac.items():
        ds = [d * p**k for d in ds for k in range(e+1)]
    return sorted(ds)

FACS = {(10,16): {13:1, 499:1}, (12,20): {5:1, 59:1, 1753:1},
        (17,27): {5:1, 71:1, 14303:1}}

for A, B in [(10, 16), (12, 20)]:
    D = (1 << B) - 3**A
    N = vectors_N(A, B, D)
    S = np.fft.ifft(N.astype(np.float64)) * D
    print("="*74)
    print(f"(A,B)=({A},{B})  D={D}")
    print("="*74)

    # Q1: orbit sums with multiplicity
    visited = np.zeros(D, dtype=bool); visited[0] = True
    sums = {}
    for t0 in range(1, D):
        if visited[t0]: continue
        t = t0; osum = 0.0+0.0j
        while not visited[t]:
            visited[t] = True; osum += S[t]; t = (2*t) % D
        key = round(osum.real, 3)
        sums.setdefault(key, []).append(t0)
    from collections import Counter
    mult = sorted(((len(v), k) for k, v in sums.items()), reverse=True)
    print("orbit-sum multiplicities (count x value), top 6:")
    for cnt, val in mult[:6]:
        print(f"   {cnt} orbits share sum {val}")
    # group-theoretic predictions
    def order(a, M):
        x = a % M; k = 1
        while x != 1: x = x*a % M; k += 1
        return k
    o2 = order(2, D)
    phiD = 1
    for p, e in FACS[(A,B)].items(): phiD *= (p-1)*p**(e-1)
    Qsize = phiD // o2
    # order of image of 3 in Q = min k with 3^k in <2>
    pow3 = 3 % D; k3 = 1
    two_powers = set()
    x = 1
    for _ in range(o2): two_powers.add(x); x = 2*x % D
    while pow3 not in two_powers:
        pow3 = pow3*3 % D; k3 += 1
    minus1 = (D-1) in two_powers
    print(f"|Q| = phi(D)/ord(2) = {phiD}/{o2} = {Qsize}; "
          f"image of 3 in Q has order {k3} (divides A={A}: {A % k3 == 0}); "
          f"-1 in <2>: {minus1}")

    # Q2: is N invariant under r -> 3r?
    idx = np.arange(D)
    inv3 = pow(3, -1, D)
    N3 = N[(inv3*idx) % D]
    print(f"N invariant under r->3r: {bool(np.array_equal(N, N3))}"
          f"   (L1 diff = {int(np.abs(N-N3).sum())} of {int(N.sum())})")

    # Q3: divisor-level counts
    Sn = comb(B-1, A-1)
    print(f"divisor-level counts Z_d = #(d | C) vs random model S/d  (S={Sn}):")
    print(f"{'d':>10} {'Z_d':>10} {'S/d':>12} {'(Z_d-S/d)/sqrt(S/d)':>20}")
    for d in divisors(D, FACS[(A,B)]):
        if d == 1: continue
        Zd = int(N[idx % d == 0].sum())
        exp = Sn/d
        z = (Zd - exp)/np.sqrt(exp) if exp > 0 else float('nan')
        print(f"{d:>10} {Zd:>10} {exp:>12.2f} {z:>20.2f}")
    print()

# (17,27): divisor counts only (heavy vector)
A, B = 17, 27
D = (1 << B) - 3**A
N = vectors_N(A, B, D)
Sn = comb(B-1, A-1)
idx = np.arange(D)
print("="*74)
print(f"(A,B)=({A},{B})  D={D}   divisor-level counts only")
print("="*74)
print(f"{'d':>10} {'Z_d':>10} {'S/d':>12} {'z-score':>10}")
for d in divisors(D, FACS[(A,B)]):
    if d == 1: continue
    Zd = int(N[idx % d == 0].sum())
    exp = Sn/d
    z = (Zd - exp)/np.sqrt(exp)
    print(f"{d:>10} {Zd:>10} {exp:>12.2f} {z:>10.2f}")
print()
print("Interpretation: |z| ~ O(1) at every divisor level = random-consistent;")
print("systematic large |z| would be exploitable arithmetic bias.")
