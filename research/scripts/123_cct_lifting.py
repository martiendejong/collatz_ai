"""
123_cct_lifting.py
==================
LIFTING STRUCTURE OF CCT ELEMENTS: How does CCT_N relate to CCT_{N+1}?

Theory: For r in CCT_N with j_old >= 1, we predict:
  - r itself maps to j_new = j_old+1 at mod-2^{N+1} (SAME r, but j shifts up)
    (only if v2(m_red*3^K - 1) = l0 is unchanged when viewed mod 2^{N+1})
  - r+2^N maps to j_new = j_old+1 at mod-2^{N+1} (LIFT, but different m_red)

This script:
1. Computes CCT_N and CCT_{N+1} explicitly
2. Checks whether each pair (K,l0) in CCT_N has a corresponding pair in CCT_{N+1}
3. Identifies the "new" CCT elements that appear only at CCT_{N+1}
4. Verifies the j-shift formula: j at mod-2^{N+1} = j at mod-2^N + 1
"""
import math

def v2(x):
    if x == 0: return 999
    c = 0
    while x % 2 == 0: x //= 2; c += 1
    return c

def modinv(a, m):
    g, x, _ = egcd(a, m)
    if g != 1: raise ValueError(f"No inverse: {a} mod {m}")
    return x % m

def egcd(a, b):
    if a == 0: return b, 0, 1
    g, x, y = egcd(b % a, a)
    return g, y - (b // a) * x, x

def cct_element(K, l0, N):
    j = N - K - l0
    if j < 1: return None
    NK = N - K; mod = 1 << NK
    target = (1 - (1 << l0)) % mod
    pow3K = pow(3, K, mod)
    inv3K = modinv(pow3K, mod)
    m_red = (target * inv3K) % mod
    if m_red % 2 == 0: return None
    r = ((1 << K) * m_red - 1) % (1 << N)
    if v2(r+1) != K: return None
    prod_exact = m_red * (3**K)
    l0_check = v2(prod_exact - 1)
    if l0_check != l0: return None
    return (r, m_red, j)

def cct_full(N):
    elements = {}  # (K,l0) -> (r, m_red, j)
    for K in range(1, N-1):
        for l0 in range(1, N-K):
            result = cct_element(K, l0, N)
            if result is not None:
                r, m_red, j = result
                elements[(K, l0)] = (r, m_red, j)
    return elements

print("=" * 70)
print("PART 1: CCT LIFTING N=8 -> N=9 -> N=10 -> N=11")
print("=" * 70)
print()

for N_src in [8, 9, 10]:
    N_dst = N_src + 1
    cct_src = cct_full(N_src)
    cct_dst = cct_full(N_dst)
    dst_by_r = {r: (K, l0, j) for (K, l0), (r, m_red, j) in cct_dst.items()}
    dst_rs = set(dst_by_r.keys())

    print(f"Lifting CCT_{N_src} (size {len(cct_src)}) -> CCT_{N_dst} (size {len(cct_dst)})")
    print(f"  New (K,l0) pairs in CCT_{N_dst}: j=1, K from 1 to N-2={N_dst-2}:")
    new_pairs = [(K, l0) for (K,l0) in cct_dst if (K,l0) not in cct_src]
    for K, l0 in sorted(new_pairs):
        r, m_red, j = cct_dst[(K, l0)]
        print(f"    (K={K}, l0={l0}, j={j}): r={r}")
    print()

    # Check how each CCT_N element maps to CCT_{N+1}
    match_count = 0
    jshift_count = 0
    print(f"  {'(K,l0)':>12} {'j_src':>6} {'r_src':>6} {'r_src in CCT_{N_dst}?':>22} {'j_dst?':>8}")
    print("  " + "-" * 60)
    for (K, l0), (r, m_red, j_src) in sorted(cct_src.items()):
        # Check if r appears in CCT_{N+1} (as an element mod 2^{N+1})
        # r mod 2^N might lift to r or r+2^N in mod-2^{N+1}
        in_dst_low = r in dst_rs
        in_dst_high = (r + (1<<N_src)) in dst_rs
        if in_dst_low:
            j_dst = dst_by_r[r][2]
            status = f"YES (r={r}, j_dst={j_dst})"
            if j_dst == j_src + 1: jshift_count += 1
        elif in_dst_high:
            j_dst = dst_by_r[r + (1<<N_src)][2]
            status = f"YES+2^N (r={r+(1<<N_src)}, j_dst={j_dst})"
            if j_dst == j_src + 1: jshift_count += 1
        else:
            status = "NOT in CCT_{N+1}"
        match_count += 1 if (in_dst_low or in_dst_high) else 0
        print(f"  ({K:>2},{l0:>2}): j={j_src:>3} r={r:>5}  {status}")
    print(f"\n  Summary: {match_count}/{len(cct_src)} CCT_{N_src} elements appear in CCT_{N_dst}")
    print(f"           {jshift_count}/{len(cct_src)} elements have j shifted by +1")
    print()

# PART 2: General lifting formula verification
print("=" * 70)
print("PART 2: LIFTING FORMULA VERIFICATION (j shifts by +1)")
print("=" * 70)
print()
print("Claim: For (K,l0) pair with N-K-l0=j>=1:")
print("  CCT element at mod-2^N has r_N")
print("  CCT element at mod-2^{N+1} for same (K,l0) has r_{N+1}")
print("  j at mod-2^{N+1} = j + 1 (since N+1-K-l0 = j+1)")
print("  r_{N+1} mod 2^N might or might NOT equal r_N")
print()
print("Showing r_N vs r_{N+1} mod 2^N:")
for N in [8, 9, 10]:
    cct_N = cct_full(N)
    cct_N1 = cct_full(N+1)
    print(f"N={N}:")
    for (K, l0), (r_N, m_red_N, j_N) in sorted(cct_N.items())[:8]:
        if (K, l0) in cct_N1:
            r_N1, m_red_N1, j_N1 = cct_N1[(K, l0)]
            match = "same" if r_N1 % (1<<N) == r_N else f"DIFF (mod 2^N: {r_N1%(1<<N)} vs {r_N})"
            print(f"  (K={K},l0={l0}): j={j_N}->{j_N1}  r={r_N}->{r_N1} mod2^N={r_N1%(1<<N)} [{match}]")
    print()

# PART 3: The NEW CCT elements at each N — their structure
print("=" * 70)
print("PART 3: STRUCTURE OF NEW CCT ELEMENTS AT EACH N")
print("=" * 70)
print()
print("New CCT elements (only at CCT_N, not at CCT_{N-1} when embedded):")
for N in [9, 10, 11, 12]:
    cct_N = cct_full(N)
    cct_prev = cct_full(N-1)
    new_pairs = [(K,l0) for (K,l0) in cct_N if (K,l0) not in cct_prev]
    print(f"N={N}: {len(new_pairs)} new pairs (all have j=1, K+l0={N-1})")
    for K, l0 in sorted(new_pairs)[:6]:
        r, m_red, j = cct_N[(K, l0)]
        print(f"  K={K:2d} l0={l0:2d} j={j} m_red={m_red:6d} r={r:6d}")
    if len(new_pairs) > 6: print(f"  ... ({len(new_pairs)-6} more)")
    print()
