# 354: cycle closure — Martien's question about the "5/8 mechanism".
# A 1-cycle (ONE climb of k ones, then j halvings, back to start) requires:
#   n = a*2^k - 1,  T^k(n) = a*3^k - 1,  (a*3^k - 1)/2^j = n
#   =>  a * (2^(k+j) - 3^k) = 2^j - 1     (the closure equation)
# So a must be (2^j - 1)/(2^(k+j) - 3^k): the "5/8-mechanism" is the demand
# that 2^(k+j) - 3^k divides 2^j - 1. Enumerate ALL solutions (pos and neg).
print("1-cycle closure equation a = (2^j-1)/(2^(k+j)-3^k), all |a| solutions k,j <= 400:")
sols = []
for k in range(1, 401):
    p3 = 3**k
    for j in range(1, 401):
        d = 2**(k + j) - p3
        if d != 0 and (2**j - 1) % d == 0:
            a = (2**j - 1) // d
            if a != 0 and a % 2 == 1:
                n = a * 2**k - 1
                sols.append((k, j, a, n))
for k, j, a, n in sols:
    print(f"  k={k}, j={j}: a={a}  ->  n = {n}   "
          f"({'TRIVIAL CYCLE 1-4-2-1' if n == 1 else 'NEGATIVE cycle' if n < 0 else 'POSITIVE NONTRIVIAL!'})")
print()

# The general cycle bound: for ANY cycle (any number of climb/descent phases m)
# with K odd steps, N halvings, minimum element n_min:
#   2^N = 3^K * prod(1 + 1/(3 n_i))  =>  0 < N ln2 - K ln3 < K/(3 n_min)
#   =>  0 < N/K - log2(3) < 1/(3 ln2 n_min)
# All correction terms are POSITIVE: extra phases cannot cancel each other;
# they only tighten the required rational approximation N/K of log2(3).
# With verification n_min = 2^68 (Barina), compute the minimal K via the
# continued fraction of log2(3).
from decimal import Decimal, getcontext
getcontext().prec = 120
log23 = Decimal(3).ln() / Decimal(2).ln()
eps = 1 / (3 * Decimal(2).ln() * Decimal(2)**68)
print(f"required one-sided approximation quality: N/K - log2(3) < {eps:.3E}")

# continued fraction convergents of log2(3)
x = log23
a0 = int(x)
cf = [a0]
frac = x - a0
h0, h1 = 1, a0
k0, k1 = 0, 1
convs = [(a0, 1)]
for _ in range(60):
    x = 1 / frac
    ai = int(x)
    frac = x - ai
    h0, h1 = h1, ai * h1 + h0
    k0, k1 = k1, ai * k1 + k0
    convs.append((h1, k1))
    if k1 > 10**14:
        break

print("convergents p/q of log2(3) with one-sided error (p/q > log2(3) only):")
best = None
for p, q in convs:
    err = Decimal(p) / Decimal(q) - log23
    side = ">" if err > 0 else "<"
    mark = ""
    if err > 0 and err < eps and best is None:
        best = (p, q)
        mark = "   <-- FIRST that satisfies the cycle requirement"
    if q > 100:
        print(f"  q={q:>15}  err={err:.3E} ({side}){mark}")
# semiconvergents (mediants) can satisfy the one-sided condition earlier than
# the next pure convergent — scan them for the sharp minimal K
cands = []
if best:
    cands.append(best[1])
for i in range(len(convs) - 1):
    p1, q1 = convs[i]
    p2, q2 = convs[i + 1]
    e1 = Decimal(p1) / Decimal(q1) - log23
    if e1 <= 0:
        continue  # need the upper-side convergent as anchor
    for t in range(1, 2000):
        p, q = p1 + t * p2, q1 + t * q2
        e = Decimal(p) / Decimal(q) - log23
        if e <= 0:
            break
        if e < eps:
            cands.append(q)
            break
K = min(cands)
length = int(Decimal(K) * (1 + log23))
print()
print(f"sharp minimal K (odd steps) for ANY cycle, incl. semiconvergents: {K:,}")
print(f"minimal total cycle length K+N ~ K*(1+log2 3) = {length:,} steps")
