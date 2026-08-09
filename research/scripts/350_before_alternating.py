# 350: what comes BEFORE the alternating family P(K) = (2^(K+2)-5)/3 = 1010...1001?
# Colour law (Obs 536): P(K) mod 3 cycles 0,2,1 over K = 3,5,7,... So:
#   K = 3 mod 6: nothing (leaf)
#   K = 5 mod 6: canonical predecessor via v=1:  Q = (2P-1)/3 = (2^(K+3)-13)/9
#   K = 1 mod 6: canonical predecessor via v=2:  Q = (4P-1)/3 = (2^(K+4)-23)/9
# Conjecture from the 2-adic frame: level-2 ancestors repeat with block length
# ord(2 mod 9) = 6, i.e. binary = (111000)^j + fixed seed tail; seeds 27 and 225.
# Level 3 should then have block length ord(2 mod 27) = 18.
def T(m):
    x = 3 * m + 1
    while x % 2 == 0:
        x //= 2
    return x

def canon_pred(n):
    if n % 3 == 0:
        return None, None
    v = 2 if n % 3 == 1 else 1
    m = (2**v * n - 1) // 3
    assert m % 2 == 1 and T(m) == n
    return m, v

print("level 2: predecessors of the alternating numbers P(K)")
for K in range(5, 42, 2):
    P = (2**(K + 2) - 5) // 3
    m, v = canon_pred(P)
    if m is None:
        print(f"  K={K:>2}: P={bin(P)[2:]}  <- niets (K=3 mod 6)")
        continue
    # closed-form check
    if K % 6 == 5:
        cf = (2**(K + 3) - 13) // 9
        fam = "A (seed 27, tail 11011)"
    else:
        cf = (2**(K + 4) - 23) // 9
        fam = "B (seed 225, tail 11100001)"
    assert m == cf, (K, m, cf)
    b = bin(m)[2:]
    # strip repeating 111000 blocks from the front
    j = 0
    while b.startswith("111000"):
        b = b[6:]
        j += 1
    print(f"  K={K:>2}: {('111000 ' * j)}{b}   family {fam}, {j} blocks, closed form OK")

print()
print("level 3: predecessors of the level-2 ancestors (expect 18-bit blocks)")
rows = []
for K in range(5, 130, 2):
    P = (2**(K + 2) - 5) // 3
    m2, _ = canon_pred(P)
    if m2 is None:
        continue
    m3, v3 = canon_pred(m2)
    if m3 is None:
        rows.append((K, m2, None, None))
    else:
        rows.append((K, m2, m3, v3))
# group by K mod 18 and check: consecutive members in a group share their tail,
# differing by one prepended 18-bit block
from collections import defaultdict
groups = defaultdict(list)
for K, m2, m3, v3 in rows:
    if m3 is not None:
        groups[K % 18].append((K, m3))
for g in sorted(groups):
    mem = groups[g]
    if len(mem) < 2:
        continue
    ok = True
    blocks = set()
    for (K1, x1), (K2, x2) in zip(mem, mem[1:]):
        b1, b2 = bin(x1)[2:], bin(x2)[2:]
        if len(b2) - len(b1) != 18 or not b2.endswith(b1):
            ok = False
            break
        blocks.add(b2[: len(b2) - len(b1)])
    K0, x0 = mem[0]
    print(f"  K mod 18 = {g:>2}: seed {x0} = {bin(x0)[2:]}")
    print(f"      18-bit prepend-block consistent: {ok}, block = {sorted(blocks) if ok else '-'}")
print()
print("leaves at level 3 (chain stops):")
for K, m2, m3, v3 in rows[:12]:
    if m3 is None:
        print(f"  K={K:>2}: level-2 ancestor {m2} = 3*{m2//3} is a leaf")
