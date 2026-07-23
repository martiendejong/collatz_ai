"""
155_affine_renewal.py
======================
THE AFFINE RENEWAL IDENTITY -- searching for individual-level structure.

The shift map sigma: (E_1,...,E_{A-1}) -> (E_1+1,...,E_{A-1}+1) is a bijection

    {patterns with E_{A-1} <= B-2}  -->  {patterns with E_1 >= 2}

and transforms C = sum_i 3^(A-1-i) 2^(E_i) EXACTLY (integer identity):

    C(sigma E) = 2*C(E) - 3^(A-1).

Define phi(r) = 2r - 3^(A-1) mod M (affine, fixed point r* = 3^(A-1)).
Then for EVERY modulus M and every residue r:

    N(phi(r)) - N_first(phi(r)) = N(r) - N_last(r)          (*)

where N(r)      = #{patterns : C == r mod M},
      N_first(r) = #{... and E_1 = 1},
      N_last(r)  = #{... and E_{A-1} = B-1}.

This is exact, proven by bijection -- no heuristics. Consequences to test:

 1. VERIFY (*) exactly at M = D = 2^B - 3^A (the real modulus).
 2. MECHANISM of Exact Uniformity mod p: iterating (*) around a phi-orbit,
    N is constant on orbits as soon as the boundary counts N_first, N_last
    are themselves uniform -- uniformity should propagate upward from
    shorter patterns. Test: find, for p=5,7, the B at which N_first/N_last
    go exactly uniform vs the B at which N does.
 3. THE ORBIT OF 0 mod D: phi^k(0) = r*(1 - 2^k). Note 2^B == 3^A (mod D),
    so the subgroup <2> mod D contains 3^A -- the affine orbit structure of
    0 is arithmetically tied to the cycle equation itself. Compute the orbit
    length of 0 under phi mod D and the telescoped expression of N(0).
"""
import numpy as np
from math import comb, log2

THETA = log2(3)

def vectors(A, B, M):
    """N, N_first (E_1=1), N_last (E_{A-1}=B-1) as exact count vectors mod M."""
    w = [pow(3, A-1-i, M) for i in range(A)]
    # dp[i][r]: i elements placed (incl E_0=0), partial C == r
    dp = [np.zeros(M, dtype=np.int64) for _ in range(A+1)]
    dp[1][w[0] % M] = 1
    # also track: dpF = restricted to E_1 = 1 (second element at position 1)
    dpF = [np.zeros(M, dtype=np.int64) for _ in range(A+1)]
    for d in range(1, B):
        p2 = pow(2, d, M)
        for i in range(min(d, A-1), 0, -1):
            c = (w[i]*p2) % M
            add = np.roll(dp[i], c)
            dp[i+1] += add
            if i >= 2:
                dpF[i+1] += np.roll(dpF[i], c)
            if i == 1 and d == 1:
                dpF[2] += np.roll(dp[1], c)      # E_1 = 1 exactly
        # last-position tracking handled separately below
    N = dp[A].copy()
    NF = dpF[A].copy()
    # N_last: E_{A-1} = B-1 -> contribution of patterns whose last chosen d = B-1.
    # Recompute dp without position B-1, then the last step must use d = B-1:
    dp2 = [np.zeros(M, dtype=np.int64) for _ in range(A+1)]
    dp2[1][w[0] % M] = 1
    for d in range(1, B-1):
        p2 = pow(2, d, M)
        for i in range(min(d, A-1), 0, -1):
            c = (w[i]*p2) % M
            dp2[i+1] += np.roll(dp2[i], c)
    cL = (w[A-1] * pow(2, B-1, M)) % M
    NL = np.roll(dp2[A-1], cL)
    return N, NF, NL

print("="*76)
print("TEST 1: EXACT VERIFICATION OF THE AFFINE RENEWAL IDENTITY MOD D")
print("="*76)
for A, B in [(5, 8), (7, 12), (10, 16), (12, 20)]:
    D = (1 << B) - 3**A
    N, NF, NL = vectors(A, B, D)
    S = comb(B-1, A-1)
    assert int(N.sum()) == S and int(NF.sum()) == comb(B-2, A-2) \
        and int(NL.sum()) == comb(B-2, A-2), "sanity totals"
    rstar = pow(3, A-1, D)
    idx = np.arange(D)
    phi = (2*idx - rstar) % D          # phi(r)
    lhs = N[phi] - NF[phi]             # N(phi(r)) - N_first(phi(r))
    rhs = N - NL                       # N(r) - N_last(r)
    ok = bool(np.array_equal(lhs, rhs))
    print(f"(A,B)=({A},{B}) D={D}:  identity holds for all {D} residues: {ok}")
print()

print("="*76)
print("TEST 2: UNIFORMITY PROPAGATION MECHANISM MOD p")
print("="*76)
print("At which B do N, N_first, N_last become EXACTLY uniform mod p?")
for p in [5, 7]:
    print(f"p = {p}:")
    for B in range(20, 301, 20):
        A = round(B/THETA)
        N, NF, NL = vectors(A, B, p)
        u  = len(set(N.tolist())) == 1
        uF = len(set(NF.tolist())) == 1
        uL = len(set(NL.tolist())) == 1
        print(f"  B={B:>3} A={A:>3}: N uniform={str(u):>5}  "
              f"N_first={str(uF):>5}  N_last={str(uL):>5}")
        if u and uF and uL and B > 100:
            break
print()

print("="*76)
print("TEST 3: THE AFFINE ORBIT OF 0 MOD D")
print("="*76)
for A, B in [(10, 16), (12, 20), (17, 27)]:
    D = (1 << B) - 3**A
    rstar = pow(3, A-1, D)
    # orbit length of 0 under phi = smallest k>0 with r*(1-2^k) == 0, i.e. 2^k == 1
    # mod D/gcd(D, r*): since gcd(r*,D)=1 (D coprime to 3), length = ord_D(2).
    k, x = 1, 2 % D
    while x != 1:
        x = 2*x % D; k += 1
        if k > 10**7: k = -1; break
    print(f"(A,B)=({A},{B}) D={D}: ord_D(2) = {k}"
          f"  (orbit of 0 under phi has this length)")
    print(f"   note: 2^B == 3^A (mod D) and B={B}, so ord_D(2) | lcm structure"
          f" linking <2>, <3> mod D")
print()
print("""INTERPRETATION
--------------
(*) is exact at the real modulus D. Telescoping it k times around the orbit:

  N(0) = N(phi^-k(0)) + sum_{j=1..k} [ N_first(phi^-j+1(0)) - N_last(phi^-j(0)) ]

Around the full orbit (k = ord_D(2)) the N-terms cancel and

  sum over the orbit of  [N_first - N_last]  =  0        (exact constraint)

The counts N_first/N_last are (A-1)-level pattern counts: the identity builds
an exact recursion tree over pattern length. NO cycle exists iff N(0) = 0 at
every admissible signature -- the recursion expresses N(0) through boundary
counts at all depths. This is genuine individual-level structure (exact, not
statistical) -- but exploiting it at depth A ~ 10^10 requires controlling the
recursion globally, which we cannot yet do.""")
