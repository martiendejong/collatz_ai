"""
159_class_invariants.py
========================
SYSTEMATIC CHARACTERIZATION OF THE <2,3>-CLASS INVARIANTS (after Obs 306).

By Orbit-Sum Invariance, D*N(0) = sum over divisor levels of a short sum of
class invariants on (Z/D)^x / <2,3>. Per signature we measure:

  - factorization of D, phi(D), ord_D(2)
  - k3 = order of the image of 3 in Q = (Z/D)^x/<2>   (must divide A, since
    3^A == 2^B mod D -- the cycle equation living inside the group theory)
  - index = [Z_D^x : <2,3>] = |Q| / k3   (number of primitive class invariants)
  - the invariant values P_class = sum of S(t) over each class, their spread,
    and the exact Ramanujan sum rule: sum_class P_class = sum_r N(r) c_D(r).

Question: is there a pattern in (k3, index, invariant values) along the
convergent line that a specialist could exploit?
"""
import numpy as np
from math import gcd, comb, log2
import random, sys

THETA = log2(3)

# ---------- factoring ----------
def is_prime(n):
    if n < 2: return False
    for p in (2,3,5,7,11,13,17,19,23,29,31,37):
        if n % p == 0: return n == p
    d, s = n-1, 0
    while d % 2 == 0: d //= 2; s += 1
    for a in (2,3,5,7,11,13,17,19,23,29,31,37):
        x = pow(a, d, n)
        if x in (1, n-1): continue
        for _ in range(s-1):
            x = x*x % n
            if x == n-1: break
        else: return False
    return True

def pollard(n):
    if n % 2 == 0: return 2
    while True:
        x = random.randrange(2, n); y = x; c = random.randrange(1, n); d = 1
        while d == 1:
            x = (x*x+c) % n; y = (y*y+c) % n; y = (y*y+c) % n
            d = gcd(abs(x-y), n)
        if d != n: return d

def factor(n):
    if n == 1: return {}
    if is_prime(n): return {n: 1}
    d = pollard(n); f = factor(d)
    for p, e in factor(n//d).items(): f[p] = f.get(p, 0)+e
    return f

def N_vector(A, B, M):
    w = [pow(3, A-1-i, M) for i in range(A)]
    dp = [np.zeros(M, dtype=np.int32) for _ in range(A+1)]
    dp[1][w[0] % M] = 1
    for d in range(1, B):
        p2 = pow(2, d, M)
        for i in range(min(d, A-1), 0, -1):
            dp[i+1] += np.roll(dp[i], (w[i]*p2) % M)
    return dp[A]

def order_mod(a, M):
    x = a % M; k = 1
    while x != 1:
        x = x*a % M; k += 1
    return k

print(f"{'(A,B)':>8} {'D':>9} {'factors':>18} {'ord2':>8} {'k3':>4} {'k3|A':>5} "
      f"{'index':>6} {'class invariants (primitive level)':>40}")
print("-"*115)

SIGS = [(5,8), (7,12), (9,15), (10,16), (11,18), (12,20), (13,21), (14,23), (15,24)]
results = []
for A, B in SIGS:
    D = (1 << B) - 3**A
    fac = factor(D)
    phi = 1
    for p, e in fac.items(): phi *= (p-1)*p**(e-1)
    o2 = order_mod(2, D)
    # k3 = order of 3 in Q: min k with 3^k in <2>
    two_pows = set(); x = 1
    for _ in range(o2): two_pows.add(x); x = 2*x % D
    y = 3 % D; k3 = 1
    while y not in two_pows:
        y = y*3 % D; k3 += 1
    index = phi // (o2 * k3)
    # class invariants
    N = N_vector(A, B, D)
    S = np.fft.ifft(N.astype(np.float64)) * D
    seen = bytearray(D)
    classes = []
    for t0 in range(1, D):
        if seen[t0] or gcd(t0, D) != 1: continue
        stack = [t0]; seen[t0] = 1; total = 0.0+0.0j; size = 0
        while stack:
            t = stack.pop()
            total += S[t]; size += 1
            for m in (2, 3):
                u = m*t % D
                if not seen[u]: seen[u] = 1; stack.append(u)
        classes.append((total.real, total.imag, size))
    # Ramanujan sum-rule check: sum over units of S = sum_r N(r) c_D(r)
    # c_D(r) = sum_{d | gcd(r,D)} d mu(D/d); compute via divisor counts
    def divisors_of(n, f):
        ds = [1]
        for p, e in f.items(): ds = [d*p**k for d in ds for k in range(e+1)]
        return ds
    def mu(n, f):
        for p, e in f.items():
            if e > 1: return 0
        return (-1)**len(f)
    idx = np.arange(D)
    ram = 0.0
    for d in divisors_of(D, fac):
        rest = D // d
        frest = factor(rest)
        m = mu(rest, frest)
        if m == 0: continue
        ram += d * m * int(N[idx % d == 0].sum())
    tot_cls = sum(c[0] for c in classes)
    ok = abs(tot_cls - ram) < 0.5
    vals = ", ".join(f"{c[0]:+.0f}" for c in sorted(classes, reverse=True)[:4])
    if len(classes) > 4: vals += ", ..."
    fs = "*".join(f"{p}" if e == 1 else f"{p}^{e}" for p, e in sorted(fac.items()))
    print(f"({A:>2},{B:>2}) {D:>9} {fs:>18} {o2:>8} {k3:>4} {str(A%k3==0):>5} "
          f"{index:>6}  [{vals}]  sumrule:{'OK' if ok else 'FAIL'}")
    results.append((A, B, D, k3, index, [c[0] for c in classes]))

print()
print("PATTERN SEARCH")
print("-"*70)
print("k3 values vs A:", [(A, k3) for A, B, D, k3, ix, _ in results])
print("index values:  ", [(f"({A},{B})", ix) for A, B, D, k3, ix, _ in results])
inv_counts = [(f"({A},{B})", len(v)) for A, B, D, k3, ix, v in results]
print("class counts:  ", inv_counts)
# normalized spread of invariants
print()
print("normalized invariants P_class / (S(0)/index) -- do they cluster?")
for A, B, D, k3, ix, vals in results:
    S0 = comb(B-1, A-1)
    if ix and vals:
        norm = [v/(S0/max(len(vals),1)) for v in sorted(vals, reverse=True)]
        print(f"  ({A:>2},{B:>2}): " + ", ".join(f"{x:+.3f}" for x in norm[:6]))
