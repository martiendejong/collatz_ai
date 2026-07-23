"""
154_equidistribution_probe.py
==============================
FIRST STEPS AT THE EQUIDISTRIBUTION WALL.

The cycle case reduces to: C = sum_i 3^(A-i) 2^(E_{i-1}) never hits 0 mod
D = 2^B - 3^A over admissible sequences (script 152). The missing ingredient
is equidistribution of C mod D. Two probes:

PART 1 -- fixed small primes p: exact residue counts N_p(r) of C mod p via
  transfer DP (exact big-integer arithmetic), along the cycle-relevant line
  A = round(B/theta). Discrepancy Delta(B) = max_r |N(r) - S/p| / S should
  decay like gamma_p^B (transfer-operator spectral gap). If it does,
  equidistribution mod any FIXED p is provable-in-principle by finite
  transfer-matrix computation -- an honest incremental result. (It cannot
  close the problem: the needed modulus D grows like 2^B.)

  Known exception p=3: C = 2^(E_{A-1}) mod 3 hits only {1,2} -- but 3 never
  divides D = 2^B - 3^A, so this costs nothing.

PART 2 -- the real modulus D: full residue vector N(r) mod D at the largest
  computable signatures. If C mod D behaved randomly, N(r) ~ Binomial(S, 1/D):
  variance/mean ~ 1, chi^2/dof ~ 1. Structure (e.g. from the rotation symmetry
  C -> 3*2^-e1*C mod D) would inflate the variance. We measure it.
"""
import numpy as np
from math import log2, comb

THETA = log2(3)

# ---------------- PART 1: fixed prime p ----------------
def residue_counts_mod_p(A, B, p):
    w = [pow(3, A-1-i, p) for i in range(A)]
    dp = [[0]*p for _ in range(A+1)]
    dp[1][w[0] % p] = 1
    for d in range(1, B):
        p2 = pow(2, d, p)
        for i in range(min(d, A-1), 0, -1):
            c = (w[i]*p2) % p
            src, dst = dp[i], dp[i+1]
            for r in range(p):
                v = src[r]
                if v:
                    dst[(r+c) % p] += v
    return dp[A]

def discrepancy(N):
    S = sum(N)
    p = len(N)
    # max_r |N_r - S/p| / S, computed exactly then floated
    worst = max(abs(N[r]*p - S) for r in range(p))
    return float((worst * 10**18) // (S*p)) / 1e18, S

print("="*76)
print("PART 1: DISCREPANCY DECAY OF C MOD p ALONG A = round(B/theta)")
print("="*76)
Bs = list(range(40, 401, 40))
for p in [3, 5, 7, 11, 13, 17, 97]:
    logs = []
    for B in Bs:
        A = round(B/THETA)
        N = residue_counts_mod_p(A, B, p)
        d, S = discrepancy(N)
        logs.append((B, log2(d) if d > 0 else None, N.count(0)))
    slope = None
    pts = [(B, l) for B, l, _ in logs if l is not None and B >= 120]
    if len(pts) >= 3:
        xs, ys = zip(*pts)
        slope = np.polyfit(xs, ys, 1)[0]
    zeros = logs[-1][2]
    tail = "  ".join(f"B={B}:{l:.1f}" if l is not None else f"B={B}:exact0"
                     for B, l, _ in logs[::3])
    gap = f"gamma_p = 2^{slope:.4f} = {2**slope:.5f}/step" if slope is not None else "n/a"
    missed = f"  [misses {zeros} residues ENTIRELY]" if zeros else ""
    print(f"p={p:>3}: log2 Delta: {tail}")
    print(f"        spectral gap fit (B>=120): {gap}{missed}")
print()
print("If all gamma_p < 1 (p != 3), C mod p equidistributes exponentially fast;")
print("equidistribution mod any fixed prime is then a finite transfer-matrix check.")
print()

# ---------------- PART 2: the real modulus D ----------------
print("="*76)
print("PART 2: RANDOMNESS TEST OF C MOD D (full residue vector)")
print("="*76)
def residue_vector_mod_D(A, B):
    D = (1 << B) - 3**A
    w = [pow(3, A-1-i, D) for i in range(A)]
    dp = [np.zeros(D, dtype=np.int64) for _ in range(A+1)]
    dp[1][w[0] % D] = 1
    for d in range(1, B):
        p2 = pow(2, d, D)
        for i in range(min(d, A-1), 0, -1):
            c = (w[i]*p2) % D
            dp[i+1] += np.roll(dp[i], c)
    return dp[A], D

print(f"{'(A,B)':>9} {'D':>10} {'S=#seq':>10} {'mean':>8} {'var/mean':>9} "
      f"{'chi2/dof':>9} {'max N':>6} {'N(0)':>5}")
for A, B in [(10, 16), (12, 20), (17, 27)]:
    N, D = residue_vector_mod_D(A, B)
    S = comb(B-1, A-1)
    assert int(N.sum()) == S
    mu = S / D
    var = float(np.var(N.astype(np.float64)))
    chi2_dof = float(((N - mu)**2 / mu).sum()) / (D - 1)
    print(f"({A:>2},{B:>2}) {D:>10} {S:>10} {mu:>8.4f} {var/mu:>9.4f} "
          f"{chi2_dof:>9.4f} {int(N.max()):>6} {int(N[0]):>5}")
print("""
Random model: var/mean ~ 1, chi2/dof ~ 1. Values near 1 mean C mod D is
statistically indistinguishable from uniform at this scale -- i.e. NO
exploitable structure; values >> 1 would reveal hidden arithmetic to use.
(The rotation symmetry C -> 3*2^-e * C mod D links residues in orbits and
is the first place any excess variance would come from.)""")
