"""
236_cycle_algebraic.py
======================
SP2A: Halvings-tree constraint satisfaction for small k.
SP2C: Mod-3 type constraints on cycle elements.

For a Collatz cycle with k odd steps, halving pattern (h_1,...,h_k):
  - h_i >= 1, sum h_i = h (determined by SP1A)
  - n0 = S / (2^h - 3^k) must be a positive odd integer
  - Each cycle element n_i must be odd and positive

Mod-3 type constraint (SP2C):
  - Odd numbers fall in three residue classes mod 3: 1, 3, 5 (mod 6) -> mod 3: 1, 0, 2
  - But odd multiples of 3 (divisible by 3) cannot be in a Collatz cycle
    (if 3|n and n is odd, then 3n+1 ≡ 1 mod 3, fine, but n = 3m means
     n is an odd multiple of 3 -> this is allowed since 3m is odd only if m is odd)
  - Actually: all positive odd integers are cycle element candidates
  - The Collatz map on odd n: T(n) = (3n+1)/2^v where v = v_2(3n+1)
  - mod-3 analysis: n ≡ 0 mod 3 => 3n+1 ≡ 1 mod 3 => T(n) ≡ 1*3^{-v} mod 3...
  - More useful: n mod 3 determines the r-type in the K-L tree:
    n ≡ 1 mod 3 => r-type 0 (in K-L: r=0 nodes, i = 3s)
    n ≡ 0 mod 3 => impossible for odd n coprime to 3... wait
    n ≡ 2 mod 3 => r-type 2 (in K-L: r=2 nodes)

  Actually in K-L indexing: odd n maps to K-L node as follows:
  - n = 2^a * m + 1 style... let me think.

  In the K-L tree at depth k, node index i corresponds to starting value
  n = i (roughly). The odd positions are i ≡ 1,5 mod 6 (= 1 mod 2).
  In K-L interleaved: r = i%3, s = i//3.
  Odd i: r=0 (i=0,3,6,...), r=1 (i=1,4,7,...), r=2 (i=2,5,8,...).
  But K-L i includes even numbers too.

  Simpler approach: look at n mod 3 for odd n:
    n ≡ 1 mod 3: e.g., 1, 7, 13, 19, 25...
    n ≡ 2 mod 3: e.g., 2 (even), 5, 11, 17, 23...
    n ≡ 0 mod 3: e.g., 3, 9, 15, 21... (odd multiples of 3)

  For Collatz: if 3|n (odd n), then 3n+1 ≡ 1 (mod 3), and v_2(3n+1) >= 1.
  So 3|n is possible in a cycle. But then gcd(n,3) = 3.

  If a cycle contains n with 3|n, then n = 3m (m odd), T(n) = (9m+1)/2^v.
  9m+1 ≡ 1 mod 3, so T(n) ≡ (1/2^v) mod 3... this gets complicated.

  Key result: any cycle element n satisfies gcd(n,3) = 1 unless it's the
  trivial structure, because:
  - If 3|n in cycle, then 3n+1 ≡ 1 mod 3
  - Next few elements can be anything mod 3
  - But the cycle equation mod 3 must be consistent

  The r-type sequence in a cycle is determined by n mod 3:
    n ≡ 1 mod 3 => K-L r=0 (since i = 3s means i ≡ 0 mod 3 ≠ 1...)

  Actually I should just track the Collatz orbit mod 3 directly.
"""
import sys
from math import log2, ceil, log

LOG2_3 = log2(3)

print("236: Cycle algebraic constraints (SP2A + SP2C)")
print("=" * 72)
sys.stdout.flush()

# ============================================================
# SP2A: Enumerate halving patterns for small k
# ============================================================
print("\n--- SP2A: Enumerate halving patterns for k=1..7 ---")
print("For each k, h = ceil(k*log2(3)), enumerate all (h_1,...,h_k) with h_i>=1, sum=h.")
print("Compute n0 = S/(2^h - 3^k). Check if n0 is an odd positive integer.")
print()

def S_from_pattern(h_vec):
    """Compute S = sum_{j=0}^{k-1} 3^j * 2^{H_j} where H_j = h_{j+1}+...+h_k."""
    k = len(h_vec)
    S = 0
    for j in range(k):
        # H_j = sum of h_{j+1},...,h_k = sum(h_vec[j+1:])
        H_j = sum(h_vec[j+1:])
        S += (3**j) * (2**H_j)
    return S

def enumerate_patterns(k, h):
    """Enumerate all (h_1,...,h_k) with h_i >= 1, sum = h."""
    if k == 1:
        yield (h,)
        return
    for h1 in range(1, h - k + 2):
        for rest in enumerate_patterns(k-1, h - h1):
            yield (h1,) + rest

total_checks = 0
cycle_candidates = []

for k in range(1, 8):
    lk = k * LOG2_3
    h = ceil(lk)
    if h == lk:
        h += 1

    gap = 2**h - 3**k
    if gap <= 0:
        print(f"k={k}: h={h}, gap = 2^{h} - 3^{k} = {gap} <= 0, skip")
        continue

    n_patterns = 0
    n_candidates = 0

    print(f"k={k}: h={h}, gap={gap}, eps={h-lk:.6f}")

    for pat in enumerate_patterns(k, h):
        n_patterns += 1
        S = S_from_pattern(pat)
        if S % gap == 0:
            n0 = S // gap
            if n0 > 0 and n0 % 2 == 1:  # positive odd integer
                # Check mod-3 type of n0
                mod3 = n0 % 3
                n_candidates += 1
                cycle_candidates.append((k, h, pat, n0, mod3))
                if n_candidates <= 5:
                    print(f"  CANDIDATE: pat={pat}, S={S}, n0={n0} (n0 mod 3 = {mod3})")
        total_checks += 1

    print(f"  Patterns checked: {n_patterns}, integer candidates: {n_candidates}")
    print()
    sys.stdout.flush()

print(f"Total patterns checked: {total_checks}")
print(f"Total cycle candidates (n0 odd integer): {len(cycle_candidates)}")
print()

if cycle_candidates:
    print("Cycle candidates found (must still verify the full orbit):")
    for k, h, pat, n0, mod3 in cycle_candidates:
        print(f"  k={k}, h={h}, pat={pat}, n0={n0}, n0 mod 3 = {mod3}")
    print()

# ============================================================
# SP2C: Mod-3 type sequence in a cycle
# ============================================================
print("--- SP2C: Mod-3 type constraint for Collatz cycles ---")
print()
print("For odd n, define r(n) = n mod 3:")
print("  r=0: n = 0 mod 3 (odd multiples of 3: 3, 9, 15, 21,...)")
print("  r=1: n = 1 mod 3 (1, 7, 13, 19, 25, 31,...)")
print("  r=2: n = 2 mod 3 (5, 11, 17, 23, 29,...)")
print()
print("Collatz step on odd n: n -> T(n) = (3n+1)/2^v2(3n+1).")
print("Mod-3 transition:")
print("  n = 0 mod 3: 3n+1 = 1 mod 3 => T(n) = 1 mod 3 (if v even) or 2 mod 3 (v odd)")
print("  n = 1 mod 3: 3n+1 = 1 mod 3 => same as above")
print("  n = 2 mod 3: 3n+1 = 1 mod 3 => same as above")
print()
print("Key: 3n+1 = 1 mod 3 always, so T(n) = (3n+1)/2^v = 2^{-v} mod 3.")
print("  v even: T(n) = 1 mod 3 (r=1)")
print("  v odd:  T(n) = 2 mod 3 (r=2)")
print()
print("So after ONE Collatz step (one halving burst h_i):")
print("  h_i even => next element = 1 mod 3 (r=1)")
print("  h_i odd  => next element = 2 mod 3 (r=2)")
print()
print("Constraint on cycle (periodic orbit):")
print("  The r-type sequence (r_0, r_1,...,r_{k-1}) is determined by (h_1,...,h_k) parity.")
print("  r_j = 1 if h_{j+1} even, else 2 (for j < k).")
print("  The cycle must start and end at the same r-type: r_0 = r_{k-1+1} = r_0.")
print("  This is automatically satisfied by the periodicity.")
print()
print("The FIRST element r_0:")
print("  r_0 is determined by h_k (last halving burst): r_0 = 1 if h_k even, else 2.")
print("  So cycles with h_k even have n0 = 1 mod 3.")
print("  Cycles with h_k odd  have n0 = 2 mod 3.")
print("  Cycles with n0 = 0 mod 3 are IMPOSSIBLE.")
print()

# Verify: for each candidate, check the pattern rule
print("Verification on SP2A candidates:")
for k, h, pat, n0, mod3 in cycle_candidates:
    h_k = pat[-1]  # last halving
    predicted_r0 = 1 if h_k % 2 == 0 else 2
    check = "OK" if predicted_r0 == mod3 or (mod3 == 0 and False) else "?"
    print(f"  n0={n0} mod3={mod3}: h_k={h_k} ({'even' if h_k%2==0 else 'odd'}) => predicted r0={predicted_r0}  {check}")

print()
print("Implication: n0 mod 3 in {1,2} only (never 0).")
print("This rules out 1/3 of candidate starting points a priori.")
print()

# Show the mod-3 orbit for the trivial cycle {1}
print("--- Trivial cycle {1}: orbit verification ---")
n = 1
print(f"  n={n}: r={n%3}, h_1=?")
# 3*1+1 = 4 = 2^2 => h=2, next = 1. Full cycle: 1 -> (3*1+1)/4 = 1.
h1 = 0
m = 3*n + 1
while m % 2 == 0:
    m //= 2
    h1 += 1
print(f"  3*{n}+1 = {3*n+1} = 2^{h1} * {m}. Next = {m}.")
print(f"  h_1 = {h1} ({'even' if h1%2==0 else 'odd'}) => r_0 of next = {1 if h1%2==0 else 2}")
print(f"  n=1: r = {1%3} = 1 (as expected for h_1={h1} even)")
print()

# Compute mod-3 orbit statistics for random k=1..8 cycle candidates with best (k,h)
print("--- SP2C: r-type distribution requirement for k-cycles ---")
print(f"{'k':>4}  {'h':>5}  {'#pat':>8}  {'n0-odd-int':>10}  {'n0≡1mod3':>10}  {'n0≡2mod3':>10}")
for k in range(1, 9):
    lk = k * LOG2_3
    h = ceil(lk)
    if h == lk:
        h += 1
    gap = 2**h - 3**k
    if gap <= 0:
        continue

    if k > 7:  # too many patterns
        print(f"  k={k}: too many patterns to enumerate (use sampling)")
        continue

    n_int = 0
    n_r1 = 0
    n_r2 = 0
    n_pat = 0
    for pat in enumerate_patterns(k, h):
        n_pat += 1
        S = S_from_pattern(pat)
        if S % gap == 0:
            n0 = S // gap
            if n0 > 0 and n0 % 2 == 1:
                n_int += 1
                if n0 % 3 == 1:
                    n_r1 += 1
                elif n0 % 3 == 2:
                    n_r2 += 1
    print(f"k={k:>3}  h={h:>5}  #pat={n_pat:>8}  int={n_int:>10}  r1={n_r1:>10}  r2={n_r2:>10}")
    sys.stdout.flush()

print()
print("done")
