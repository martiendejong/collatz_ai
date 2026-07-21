"""
128_vadic_3K.py
================
The 2-adic valuation identity: v2(3^K - 1) = v2(K) + 2 for even K; = 1 for odd K.

This identity (a consequence of LTE and the order of 3 in Z/2^n) has deep
implications:
1. For n = 2^K - 1 (all-1s binary, the "extremal" Collatz element), the
   macro-step gives l0 = v2(3^K - 1) which follows this pattern.
2. This determines the DIAGONAL STRUCTURE of the CCT elements at K even vs odd.
3. It explains why K odd vs even elements in BSet behave differently.

Key questions:
1. Verify the identity numerically for K=1..20.
2. Derive the macro-step chain for n=2^K-1 (the all-1s starting points).
3. Find the ORBIT structure: 2^K-1 -> 2^{K'}-1 -> ... does this define a
   specific sub-orbit?
4. What is the period-2 structure in the BSet A/B partition in terms of v2(3^K-1)?
5. Does the identity v2(3^K-1) = v2(K)+2 connect to the CCT formula
   |CCT_N(j>=1)| = (N-2)(N-1)/2?
"""
import math
from collections import defaultdict
import numpy as np

def v2(x):
    if x == 0: return 999
    c = 0
    while x % 2 == 0: x //= 2; c += 1
    return c

def macro_step(n):
    K = v2(n+1); m = (n+1) >> K; x = m*(3**K)-1; l = v2(x)
    return x >> l, K, l

print("=" * 70)
print("PART 1: NUMERICAL VERIFICATION OF v2(3^K - 1) = v2(K)+2 (even) / 1 (odd)")
print("=" * 70)
print()
print(f"{'K':>4} {'3^K-1':>20} {'v2(3^K-1)':>12} {'v2(K)+2':>10} {'formula':>10}")
print("-" * 60)
for K in range(1, 25):
    val = 3**K - 1
    v = v2(val)
    formula = v2(K) + 2 if K % 2 == 0 else 1
    match = "OK" if v == formula else f"MISMATCH"
    print(f"{K:>4} {val:>20d} {v:>12d} {formula:>10} {match:>10}")

print()
print("=" * 70)
print("PART 2: ORBIT OF n = 2^K - 1 (THE ALL-1s ELEMENT)")
print("=" * 70)
print()
print("Starting from n=2^K-1 (K trailing 1-bits in n+1=2^K, m=1):")
print("macro_step: x = 1 * 3^K - 1 = 3^K - 1,  l0 = v2(3^K - 1)")
print("n_out = (3^K - 1) / 2^{l0}")
print()
print("Orbit chain: 2^K-1 -> n_out -> ... (track K values)")
print()

for K_start in [1, 2, 3, 4, 5, 6, 7, 8]:
    n = 2**K_start - 1
    print(f"K={K_start}: n = 2^{K_start}-1 = {n}")
    orbit = [n]
    K_orbit = [K_start]
    for step in range(12):
        n_out, K_val, l_val = macro_step(n)
        K_next = v2(n_out + 1)
        print(f"  step {step+1}: K={K_val}, l0={l_val}, n_out={n_out}, K_next={K_next}", end="")
        if n_out == 1:
            print(f" -> 1 (done)")
            break
        # Check if n_out+1 is a power of 2 (all-1s again)
        if (n_out + 1) & n_out == 0:  # power of 2
            print(f" [n_out=2^{K_next}-1, ALL-1s!]")
        else:
            print()
        n = n_out
        K_orbit.append(K_next)
        orbit.append(n)
        if step >= 11:
            print(f"  ...")
    print()

print("=" * 70)
print("PART 3: THE FRACTIONAL CASCADE — HOW DOES 3^K/2^{v2(3^K-1)} LOOK?")
print("=" * 70)
print()
print("After one macro-step from n=2^K-1 (m=1): n_out = (3^K - 1) / 2^{l0}")
print("Let's write n_out = (3^K - 1) / 2^{l0} and find v2(n_out + 1)")
print()
print(f"{'K':>4} {'l0':>4} {'n_out':>12} {'K_next':>8} {'n_out+1':>12} {'factor':>20}")
for K in range(1, 20):
    l0 = v2(3**K - 1)
    n_out = (3**K - 1) >> l0
    K_next = v2(n_out + 1)
    n_out_plus1 = n_out + 1
    # Factor n_out+1 = 2^K_next * m_next
    m_next = n_out_plus1 >> K_next
    print(f"{K:>4} {l0:>4} {n_out:>12d} {K_next:>8} {n_out_plus1:>12d} = 2^{K_next} * {m_next}")

print()
print("=" * 70)
print("PART 4: THE (K, l0) PAIR PATTERN IN THE ALL-1s ORBIT")
print("=" * 70)
print()
print("Pattern of (K, l0) pairs along orbit of n=2^30-1:")
n = 2**30 - 1
K_l0_pairs = []
for _ in range(50):
    K = v2(n+1)
    m_val = (n+1) >> K
    x = m_val * (3**K) - 1
    l0 = v2(x)
    K_l0_pairs.append((K, l0))
    n = x >> l0
    if n < 2: break

# Show the sequence of (K, l0) pairs
print("(K, l0) pairs along orbit of 2^30-1:")
for i, (K, l0) in enumerate(K_l0_pairs):
    lyap = K*math.log(3) - (K+l0)*math.log(2)
    print(f"  step {i+1:2d}: K={K:2d}, l0={l0:2d}, lyap={lyap:+.3f}")

print()
print("=" * 70)
print("PART 5: CONNECTING v2(3^K-1) TO CCT STRUCTURE")
print("=" * 70)
print()
print("CCT element with K=K, l0=l0 at modulus 2^N exists iff N > K+l0.")
print("For the 'base case' m=1 CCT element: m_red=1, so K*log3 - l0*log2 is the Lyapunov.")
print("The BASE CCT element has l0 = v2(m_red * 3^K - 1) = v2(3^K - 1) for m=1.")
print("So l0_base(K) = 1 if K odd, v2(K)+2 if K even.")
print()
print("Base CCT element (m_red=1) birth generation N0 = K + l0_base + 1:")
print(f"{'K':>4} {'l0_base':>8} {'N0=K+l0+1':>10} {'Lyapunov':>12}")
for K in range(1, 16):
    l0_base = 1 if K % 2 == 1 else v2(K) + 2
    N0 = K + l0_base + 1
    lyap = K*math.log(3) - (K+l0_base)*math.log(2)
    print(f"{K:>4} {l0_base:>8} {N0:>10d} {lyap:>12.4f}")

print()
print("Observation: for K even, l0_base = v2(K)+2 grows with K.")
print("For K=2: l0=3, N0=6. For K=4: l0=4, N0=9. For K=8: l0=5, N0=14.")
print("The base CCT element is BORN at a LATER generation for even K.")
print("This creates an ASYMMETRY between even and odd K in the CCT hierarchy.")

print()
print("=" * 70)
print("PART 6: ODD K vs EVEN K IN BSET")
print("=" * 70)
BSet = [27,55,63,83,95,103,127,159,169,191,207,223,239,253,255]
print("BSet elements by K parity:")
odd_K = [(r, v2(r+1)) for r in BSet if v2(r+1) % 2 == 1]
even_K = [(r, v2(r+1)) for r in BSet if v2(r+1) % 2 == 0]
print(f"Odd K: {odd_K}")
print(f"Even K: {even_K}")
print()

# For each BSet element, compute l0 and compare to l0_base
print("Comparing observed E[l0] with l0_base for BSet elements:")
print(f"{'r':>5} {'K':>4} {'l0_base':>8} {'E[l0]obs':>10} {'same?':>8}")
# We need the l0 data from script 125. Let's recompute here.
MOD = 256
for r in sorted(BSet):
    K = v2(r+1)
    l0_base = 1 if K % 2 == 1 else v2(K) + 2
    # Quick estimate of E[l0] by sampling
    l0_sum = 0; cnt = 0
    for k in range(500):
        n = r + MOD * k
        if v2(n+1) != K: continue
        m = (n+1) >> K; x = m*(3**K)-1; l0_sum += v2(x); cnt += 1
    E_l0 = l0_sum/cnt if cnt > 0 else 0
    same = "YES" if abs(E_l0 - l0_base) < 0.3 else "no"
    print(f"{r:>5} {K:>4} {l0_base:>8} {E_l0:>10.3f} {same:>8}")
