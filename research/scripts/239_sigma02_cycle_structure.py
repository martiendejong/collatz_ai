"""
239_sigma02_cycle_structure.py
==============================
Cycle structure of sigma_0 and sigma_2 (compared to sigma_1 which is a single Nl-cycle).

The five K-L index maps at scale Nl = 3^{k-2}:
  sigma_0(s) = 4s mod Nl          (T4-pullback for r=0)
  sigma_1(s) = (4s+2) mod Nl      (T4-pullback for r=1) -- SINGLE Nl-CYCLE (proved)
  sigma_2(s) = (4s+3) mod Nl      (T4-pullback for r=2)... wait let me re-derive.

Actually: T4(i) = (4i+2) mod N. For i = 3s+r:
  T4(3s+r) = (12s + 4r + 2) mod N = 3*((4s + floor((4r+2)/3)) mod Nl) + ((4r+2) mod 3)
For r=0: 4r+2=2, floor=0, mod3=2: T4(3s+0) = 3*(4s mod Nl) + 2. So T4 maps r=0 to r=2 with sigma_0(s)=4s.
For r=1: 4r+2=6, floor=2, mod3=0: T4(3s+1) = 3*((4s+2) mod Nl) + 0. So T4 maps r=1 to r=0 with sigma_1(s)=(4s+2).
For r=2: 4r+2=10, floor=3, mod3=1: T4(3s+2) = 3*((4s+3) mod Nl) + 1. So T4 maps r=2 to r=1 with sigma_2(s)=(4s+3).

So the chain is: r=0 --sigma_0--> r=2 --sigma_2--> r=1 --sigma_1--> r=0 (as sub-maps from s to s')
Full orbit: sigma_total = sigma_1 o sigma_2 o sigma_0 (composition on s-coords)
= sigma_1(sigma_2(sigma_0(s))) = sigma_1(sigma_2(4s)) = sigma_1(4(4s)+3) = sigma_1(16s+3) = 4(16s+3)+2 = 64s+14

sigma_total(s) = (64s+14) mod Nl.

Questions:
1. Is sigma_0 a single cycle? sigma_0(s) = 4s mod Nl. This is multiplication by 4 mod 3^{k-2}.
   Order of 4 mod 3^{k-2}: ord_{3^{k-2}}(4) = 3^{k-3}/2 if... actually 4=2^2, ord(2) mod 3^j = 2*3^{j-1}.
   So ord(4) mod 3^{k-2} = 3^{k-3} (half of ord(2)). This is NOT Nl = 3^{k-2} in general.
   => sigma_0 is NOT a single cycle (has multiple cycles of length 3^{k-3}).

2. Is sigma_2 a single cycle? sigma_2(s) = (4s+3) mod Nl.
   This is an affine map. If gcd(4, Nl) = 1 (yes, since Nl = 3^{k-2} is odd and gcd(4,3)=1),
   then sigma_2 is a permutation. But is it a single cycle?
   sigma_2^n(s) = (4^n * s + 3*(4^{n-1}+...+1)) mod Nl = 4^n*s + 3*(4^n-1)/3 mod Nl
   = 4^n*s + (4^n-1) mod Nl = (4^n+1)*s + ... wait:
   sigma_2(s) = 4s+3. sigma_2^2(s) = 4(4s+3)+3 = 16s+15. sigma_2^3(s) = 4(16s+15)+3 = 64s+63.
   sigma_2^n(s) = 4^n*s + 3*(4^{n-1}+...+1) = 4^n*s + (4^n-1) mod Nl.
   Fixed point: sigma_2(s*)=s* => 4s*+3=s* => 3s*=-3 => s*=-1 mod Nl = Nl-1.
   sigma_2 is affine with fixed point at Nl-1. Not necessarily a single cycle.

3. sigma_total(s) = (64s+14) mod Nl.
   Fixed point: 64s+14=s => 63s=-14 => s=-14/63 = -2/9 mod Nl.
   Nl = 3^{k-2}. gcd(63, Nl): 63=9*7, gcd(9*7, 3^{k-2}) = 9 for k>=4. So 63*s=-14 has solutions mod Nl/9.
   => sigma_total has multiple orbits.

This script measures the cycle structure of sigma_0, sigma_2, sigma_total.
"""
import sys
import numpy as np
from math import log2

ALPHA = log2(3.0)
LAM = 1.70
N_ITER = 500

print("239: Cycle structure of sigma_0, sigma_2, sigma_total (vs sigma_1 = single cycle)")
print(f"     lambda={LAM}")
print("=" * 72)
sys.stdout.flush()

def count_cycles(perm):
    """Count cycles and their lengths in a permutation given as array."""
    n = len(perm)
    visited = np.zeros(n, dtype=bool)
    cycle_lengths = []
    for start in range(n):
        if not visited[start]:
            cycle_len = 0
            s = start
            while not visited[s]:
                visited[s] = True
                s = perm[s]
                cycle_len += 1
            cycle_lengths.append(cycle_len)
    return sorted(cycle_lengths, reverse=True)

for k in range(4, 15):
    N  = 3 ** (k - 1)
    Nl = N // 3

    sl = np.arange(Nl, dtype=np.int64)
    sigma0 = (4 * sl) % Nl
    sigma1 = (4 * sl + 2) % Nl
    sigma2 = (4 * sl + 3) % Nl
    sigma_total = (64 * sl + 14) % Nl

    # sigma_20 (used in ve2 analysis)
    sigma_20 = (16 * sl + 14) % Nl

    cycles0 = count_cycles(sigma0)
    cycles1 = count_cycles(sigma1)
    cycles2 = count_cycles(sigma2)
    cycles_total = count_cycles(sigma_total)
    cycles_20 = count_cycles(sigma_20)

    n0 = len(cycles0); l0_max = max(cycles0) if cycles0 else 0
    n1 = len(cycles1); l1_max = max(cycles1) if cycles1 else 0
    n2 = len(cycles2); l2_max = max(cycles2) if cycles2 else 0
    nt = len(cycles_total); lt_max = max(cycles_total) if cycles_total else 0
    n20 = len(cycles_20); l20_max = max(cycles_20) if cycles_20 else 0

    print(f"k={k:2d} Nl={Nl:7d}: "
          f"sig0: {n0:4d} cycles (max len {l0_max:7d}),  "
          f"sig1: {n1:4d} cycle (len {l1_max:7d}),  "
          f"sig2: {n2:4d} cycles (max len {l2_max:7d}),  "
          f"sig_tot: {nt:4d} cyc (max {lt_max:7d}),  "
          f"sig_20: {n20:4d} cyc (max {l20_max:7d})")
    sys.stdout.flush()

print()
print("Legend:")
print("  sig1: T4-pullback for r=1 nodes, s -> (4s+2) mod Nl. PROVED single cycle (lem:sigma1).")
print("  sig0: T4-pullback for r=0, s -> 4s mod Nl. Multiplication, multiple orbits.")
print("  sig2: T4-pullback for r=2, s -> (4s+3) mod Nl. Affine, fixed pt at Nl-1.")
print("  sig_total = sig1 o sig2 o sig0: full T4^3 orbit on s-coords.")
print("  sig_20 = composed map used in ve2 analysis (sigma_20 = 16s+14).")
print()

# Verify sigma1 single-cycle claim via LTE
print("sigma1 single-cycle verification via LTE:")
print("  ord_{3^{k-2}}(4) where 4 ≡ 1+3 (one 3-unit step above 1):")
print("  LTE: v_3(4^{3^{k-2}} - 1) = v_3(4-1) + v_3(3^{k-2}) = 1 + (k-2) = k-1")
print("  sigma1(s) = 4s+2; sigma1^{Nl}(0) = 4^{Nl}*0 + 2*(4^{Nl}-1)/3 mod Nl")
print("  At s=0: sigma1^n(0) = 2*(4^n-1)/3 mod Nl.")
print("  First return to 0: 2*(4^n-1)/3 = 0 mod Nl => 4^n = 1 mod Nl = 3^{k-2}.")
print("  ord(4) mod 3^{k-2} = 3^{k-3} (half of ord(2)=2*3^{k-3}).")
print("  WAIT: this gives 3^{k-3} != Nl = 3^{k-2}. The single-cycle claim needs re-check.")
print()

# Direct check: sigma_1 orbit length starting from 0
for k in [5, 7, 9]:
    Nl = 3**(k-2)
    sl_arr = np.arange(Nl, dtype=np.int64)
    sigma1 = (4*sl_arr+2) % Nl
    # Orbit of 0 under sigma1
    s = 0
    length = 0
    visited_set = set()
    while s not in visited_set:
        visited_set.add(s)
        s = int(sigma1[s])
        length += 1
    print(f"  k={k}: Nl={Nl}, orbit of 0 under sigma1 has length {length} (Nl={Nl})")
    sys.stdout.flush()

print()
print("done")
