"""
156_renewal_consequences.py
============================
CONSEQUENCES OF THE AFFINE RENEWAL IDENTITY (continuing script 155).

 P1 (Fixed-Point Corollary). phi(r) = 2r - 3^(A-1) mod D fixes r* = 3^(A-1).
    The renewal identity at r = r* forces, EXACTLY:
        N_first(r*) = N_last(r*).
    (Patterns with E_1 = 1 hitting C == 3^(A-1) are equinumerous with patterns
    with E_{A-1} = B-1 hitting the same residue.) Verify at several signatures.

 P2 (Rotation-Divisibility Dichotomy). Rotation acts on solution patterns:
    C_rot = (3C + D)/2^(e_1) exactly, so C == 0 => C_rot == 0. A primitive
    cycle contributes exactly A distinct patterns (its rotations). Hence
        N(0) = sum_{d | gcd(A,B)} (A/d) * P(A/d, B/d),
    P = number of primitive cycles. At any signature: N(0) = 0 or
    N(0) >= A/gcd(A,B). At the minimal admissible signature A ~ 9e9:
    the number of solution patterns is ZERO or >= ~10^9. Verify the formula
    on the trivial-multiple signatures (k, 2k).

 P3 (Orbit structure). ord_D(2), ord_D(3) for near-convergent D, and the
    forced relation from 2^B == 3^A (mod D): 3^A lies in <2>, so
    ord(3)/gcd(ord(3), A) divides ord(2). Equidistribution + autocorrelation
    of N along the phi-orbit of 0 (any exploitable structure would show up
    as autocorrelation; random placement would show none).
"""
import numpy as np
from math import gcd, comb, log2
import random

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
            c = (w[i]*p2) % M
            dp2[i+1] += np.roll(dp2[i], c)
    NL = np.roll(dp2[A-1], (w[A-1]*pow(2, B-1, M)) % M)
    return dp[A], dpF[A], NL

print("="*76)
print("P1: FIXED-POINT COROLLARY  N_first(r*) = N_last(r*),  r* = 3^(A-1) mod D")
print("="*76)
for A, B in [(5, 8), (7, 12), (10, 16), (12, 20)]:
    D = (1 << B) - 3**A
    N, NF, NL = vectors(A, B, D)
    rstar = pow(3, A-1, D)
    print(f"(A,B)=({A},{B}) D={D}: N_first(r*)={int(NF[rstar])} "
          f"N_last(r*)={int(NL[rstar])}  equal: {int(NF[rstar])==int(NL[rstar])}"
          f"   [N(r*)={int(N[rstar])}]")
print()

print("="*76)
print("P2: ROTATION-DIVISIBILITY  N(0) = sum_(d|gcd) (A/d) P(A/d,B/d)")
print("="*76)
def N0(A, B):
    D = (1 << B) - 3**A
    w = [pow(3, A-1-i, D) for i in range(A)]
    dp = [np.zeros(D, dtype=np.int64) for _ in range(A+1)]
    dp[1][w[0] % D] = 1
    for d in range(1, B):
        p2 = pow(2, d, D)
        for i in range(min(d, A-1), 0, -1):
            dp[i+1] += np.roll(dp[i], (w[i]*p2) % D)
    return int(dp[A][0])

# trivial-multiple check: (k,2k) has exactly the d=k repetition of (1,2)
for k in [1, 2, 3, 4, 5]:
    A, B = k, 2*k
    n0 = N0(A, B)
    # formula: sum over d | gcd(k,2k)=k of (k/d)*P(k/d, 2k/d); P(j,2j)=1 iff j==1 else 0
    pred = 1  # only d=k contributes: (k/k)*P(1,2) = 1
    print(f"(A,B)=({k},{2*k}): N(0)={n0}, formula predicts {pred}, match={n0==pred}")
g = gcd(8963457697, 14206744327)
print(f"\nminimal admissible signature: gcd(A_min, B_min) = {g}")
print(f"=> at the minimal signature, N(0) = 0 or N(0) >= {8963457697//g:,}")
print("   (the solution count is all-or-nothing: zero, or ~10^9 patterns)")
print()

print("="*76)
print("P3: ORBIT STRUCTURE AND AUTOCORRELATION ALONG THE ORBIT OF 0")
print("="*76)
def order(a, M, cap=10**8):
    x = a % M; k = 1
    while x != 1:
        x = x*a % M; k += 1
        if k > cap: return None
    return k

for A, B in [(10, 16), (12, 20)]:
    D = (1 << B) - 3**A
    o2, o3 = order(2, D), order(3, D)
    rel = o3 // gcd(o3, A)
    print(f"(A,B)=({A},{B}) D={D}: ord(2)={o2} ord(3)={o3} "
          f"ord(3)/gcd(ord3,A)={rel} divides ord(2): {o2 % rel == 0}")

# autocorrelation of N along the phi-orbit of 0 at (12,20)
A, B = 12, 20
D = (1 << B) - 3**A
N, NF, NL = vectors(A, B, D)
rstar = pow(3, A-1, D)
orbit = []
r = 0
for _ in range(order(2, D)):
    orbit.append(r)
    r = (2*r - rstar) % D
vals = N[orbit].astype(np.float64)
t = len(vals)
mu, sd = vals.mean(), vals.std()
print(f"\n(12,20): phi-orbit of 0 has length {t}; N along orbit: "
      f"mean={mu:.4f} (global {N.mean():.4f}), sd={sd:.4f}")
# equidistribution of the orbit itself over [0, D)
srt = np.sort(np.array(orbit, dtype=np.float64))/D
ks = float(np.max(np.abs(srt - (np.arange(1, t+1))/t)))
print(f"orbit equidistribution: KS statistic = {ks:.5f}  "
      f"(random-uniform expectation ~ {1.22/np.sqrt(t):.5f})")
# autocorrelation at small lags vs shuffled null
ac = [float(np.corrcoef(vals[:-k], vals[k:])[0,1]) for k in (1, 2, 5, 10, 50)]
rng = np.random.default_rng(1)
sh = vals.copy(); rng.shuffle(sh)
acs = [float(np.corrcoef(sh[:-k], sh[k:])[0,1]) for k in (1, 2, 5, 10, 50)]
print(f"autocorr N(r_k) lags (1,2,5,10,50):  orbit   {['%+.4f'%a for a in ac]}")
print(f"                                     shuffled{['%+.4f'%a for a in acs]}")
print(f"null threshold ~ 2/sqrt(t) = {2/np.sqrt(t):.4f}")
