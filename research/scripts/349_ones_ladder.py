# 349: the predecessor ladder of the all-ones numbers 2^K - 1 (Syracuse odd-tree).
# Backward Syracuse: odd n has odd predecessors m_v = (2^v * n - 1)/3 for every v
# with 2^v * n = 1 mod 3, i.e.:
#   n = 0 mod 3  ->  NO predecessor (leaf, "niets")
#   n = 1 mod 3  ->  v even (canonical v=2)
#   n = 2 mod 3  ->  v odd  (canonical v=1)
# Martien's list: 1 -> 100 -> 1;  niets -> 11;  1001 -> 111;  niets -> 1111;  ? -> 11111.
def T(m):  # Syracuse
    x = 3 * m + 1
    while x % 2 == 0:
        x //= 2
    return x

def canon_pred(n):
    if n % 3 == 0:
        return None, None
    v = 2 if n % 3 == 1 else 1
    m = (2**v * n - 1) // 3
    assert (2**v * n - 1) % 3 == 0 and m % 2 == 1 and T(m) == n
    return m, v

print("the all-ones spine 2^K - 1:")
for K in range(1, 16):
    n = 2**K - 1
    m, v = canon_pred(n)
    if K == 1:
        print(f"  K={K:>2}  {bin(n)[2:]:>16}  (trivial cycle 1 -> 100 -> 1)")
        continue
    if m is None:
        print(f"  K={K:>2}  {bin(n)[2:]:>16}  <- niets   (K even => 2^K-1 = 0 mod 3)")
    else:
        print(f"  K={K:>2}  {bin(n)[2:]:>16}  <- {m} = {bin(m)[2:]}   (v={v}; "
              f"closed form (2^(K+2)-5)/3 = {(2**(K+2)-5)//3})")
print()
print("canonical chains backwards until a leaf (m = 0 mod 3):")
for K in [3, 5, 7, 9, 11, 13, 15]:
    n = 2**K - 1
    chain = [n]
    while True:
        m, v = canon_pred(chain[-1])
        if m is None:
            break
        chain.append(m)
        if len(chain) > 30:
            break
    s = " <- ".join(f"{x}" for x in chain)
    print(f"  K={K:>2}: {s}   (leaf: {chain[-1]} = 3*{chain[-1]//3})")
print()
print("mod-3 colour of the canonical predecessor P(K) = (2^(K+2)-5)/3, K odd:")
for K in range(3, 26, 2):
    P = (2**(K+2) - 5) // 3
    c = P % 3
    tag = {0: "leaf (niets ervoor)", 1: "continues via v even", 2: "continues via v odd"}[c]
    print(f"  K={K:>2}: P = {P:>10} = {bin(P)[2:]:>26}  P mod 3 = {c}  -> {tag}")
