# 351: can the block hierarchy DETERMINE the longest orbits?
# Test 1: are the actual total-stopping-time record holders up to 2e6 members of
#         the ancestor tree (trailing ones-layer + 10-blocks / 111000-blocks)?
# Test 2: constructor duel per bit-size: all-ones 2^B-1 vs level-2 family-A member
#         vs the true record in that dyadic range.
import sys

def steps(n, cache):
    # total stopping time (standard Collatz steps to reach 1), small-value cache
    s = 0
    m = n
    while m != 1:
        if m < len(cache) and cache[m] >= 0:
            return s + cache[m]
        m = 3 * m + 1 if m % 2 else m // 2
        s += 1
    return s

N = 2 * 10**6
cache = [-1] * (N + 1)
cache[1] = 0
for i in range(2, N + 1):
    # iterative fill
    path = []
    m = i
    while m != 1 and not (m <= N and cache[m] >= 0):
        path.append(m)
        m = 3 * m + 1 if m % 2 else m // 2
    base = cache[m] if m <= N else steps(m, cache)
    for j, x in enumerate(reversed(path)):
        if x <= N:
            cache[x] = base + j + 1

records = []
best = -1
for n in range(2, N + 1):
    if cache[n] > best:
        best = cache[n]
        records.append((n, best))

def decompose(n):
    # strip trailing ones, then greedily read 10-blocks and 111000-blocks upward
    b = bin(n)[2:]
    k = len(b) - len(b.rstrip('1'))
    rest = b[: len(b) - k] if k else b
    alt = 0
    r = rest
    # after the ones a zero must come; check for alternating (01)^m above
    while r.endswith('01'):
        r = r[:-2]
        alt += 1
    blocks6 = 0
    r2 = rest
    while r2.endswith('111000') or r2.endswith('000111'):
        r2 = r2[:-6]
        blocks6 += 1
    return b, k, alt, blocks6

print(f"total-stopping-time records up to {N} and their tree signature:")
print(f"{'n':>9} {'steps':>6}  {'binary':>24}  ones-run  10-blocks  111000-blocks")
for n, s in records[-14:]:
    b, k, alt, b6 = decompose(n)
    print(f"{n:>9} {s:>6}  {b:>24}  {k:>8}  {alt:>9}  {b6:>13}")
print()

# significance: expected ones-run of a random odd number is ~2 (geometric);
# expected 10-block count above it ~1. Sum k+2*alt is the 'tree depth in bits'.
import random
random.seed(1)
rand_depth = []
for _ in range(20000):
    n = random.randrange(3, N, 2)
    _, k, alt, _ = decompose(n)
    rand_depth.append(k + 2 * alt)
rec_depth = []
for n, s in records[-14:]:
    _, k, alt, _ = decompose(n)
    rec_depth.append(k + 2 * alt)
print(f"tree-tail depth (ones + 2*altblocks): records mean {sum(rec_depth)/len(rec_depth):.1f} bits "
      f"vs random odd mean {sum(rand_depth)/len(rand_depth):.1f} bits")
print()

# Test 2: constructor duel around B = 17 bits
def full_steps(n):
    s = 0
    while n != 1:
        n = 3 * n + 1 if n % 2 else n // 2
        s += 1
    return s

allones = 2**17 - 1
famA = int('11100011100011011', 2)  # level-2 family A, j=2 (17 bits)
# true record in [2^16, 2^17)
best_n, best_s = 0, -1
for n in range(2**16, 2**17):
    if cache[n] > best_s:
        best_s = cache[n]
        best_n = n
print("constructor duel at 17 bits:")
print(f"  all-ones 2^17-1 = {allones}: {full_steps(allones)} steps")
print(f"  level-2 family-A member {famA} = {bin(famA)[2:]}: {full_steps(famA)} steps")
print(f"  true record in [2^16,2^17): {best_n} = {bin(best_n)[2:]}: {best_s} steps")
bb, kk, aa, b66 = decompose(best_n)
print(f"    its signature: ones-run {kk}, 10-blocks {aa}, 111000-blocks {b66}")
