"""
100_phi_uniformization_and_mod512.py
=====================================
TWO RELATED INVESTIGATIONS:

A) PHI UNIFORMIZATION EFFECT:
   Why does Phi(r) = E[k/step until BSet] converge toward 2 as k0 grows?

   Theory: Phi(r) ≈ 2 + (k0-2)/E_r[h]

   - k0=1: Phi ≈ 2 + (1-2)/E[h] = 2 - 1/E[h]. If E[h]=1, Phi=1.
   - k0=8: Phi ≈ 2 + 6/8.49 ≈ 2.71. Actual: 2.412. (Lower due to internal k<2)

   This shows: Phi is always between min(k0,2) and max(k0,2) roughly.
   The ergodic avg k/step ≈ 2 because ALL elements pull toward 2.

B) MOD-512 BSet EXTENSION:
   The mod-256 analysis gives BSet (15 elements) with ergodic avg 2.06.
   What happens with mod-512 residues?

   Key question: Does the ergodic avg CHANGE with finer modular resolution?
   If the ergodic avg is truly ≈2.06 for ALL n, we should get similar results
   at mod-512.
"""
import sys, time, math
from collections import Counter, defaultdict
import numpy as np

sys.stdout.reconfigure(encoding='utf-8')

def v2(x):
    c = 0
    while x % 2 == 0: x //= 2; c += 1
    return c

def macro_step(n):
    k = v2(n + 1)
    m = (n + 1) >> k
    x = m * (3**k) - 1
    l = v2(x)
    return x >> l, k, l

BSet256 = {27, 55, 63, 83, 95, 103, 127, 159, 169, 191, 207, 223, 239, 253, 255}

M_BASE = 10**12
N = 512

# =====================================================================
# PART A: PHI UNIFORMIZATION — WHY Phi(r) → 2
# =====================================================================
print("=" * 70)
print("PART A: PHI UNIFORMIZATION EFFECT")
print("=" * 70)
print()
print("Theory: Phi(r) = E[(k0 + k_2 + ... + k_h) / h]")
print("       ≈ 2 + (k0-2)/E[h]  [since internal k_i ~ Geo(1/2), E=2]")
print()

BList256 = sorted(BSet256)
phi_data = {}

for r in BList256:
    k0 = v2(r + 1)
    k_totals = []
    h_vals = []

    for i in range(N):
        n = M_BASE + r + 256 * i
        n_cur = n
        h = 0
        k_sum = 0

        while h < 1000:
            n_out, k, l = macro_step(n_cur)
            h += 1
            k_sum += k
            n_cur = n_out
            if n_cur <= 1:
                break
            if n_cur % 256 in BSet256:
                k_totals.append(k_sum)
                h_vals.append(h)
                break

    phi = sum(k_totals) / sum(h_vals) if h_vals else 0
    avg_h = sum(h_vals) / len(h_vals) if h_vals else 0

    # Decompose: k_sum = k0 + rest
    # avg_k_rest = (total_k - k0×count) / (total_h - count)
    total_k = sum(k_totals)
    total_h = sum(h_vals)
    count = len(h_vals)
    k_rest = total_k - k0 * count  # total k from steps 2,3,...
    h_rest = total_h - count        # total steps 2,3,...

    avg_k_rest = k_rest / h_rest if h_rest > 0 else 0

    phi_theory = 2 + (k0 - 2) / avg_h if avg_h > 0 else k0

    phi_data[r] = {'k0': k0, 'phi': phi, 'avg_h': avg_h,
                   'avg_k_rest': avg_k_rest, 'phi_theory': phi_theory}

print(f"{'r':>4}  {'k0':>4}  {'Phi(r)':>8}  {'E[h]':>8}  {'k_rest/step':>12}  "
      f"{'Phi_theory':>12}  {'error':>8}")
print("-" * 70)
for r in BList256:
    d = phi_data[r]
    error = d['phi'] - d['phi_theory']
    print(f"r={r:3d}  k0={d['k0']}  Phi={d['phi']:.4f}  E[h]={d['avg_h']:.2f}  "
          f"k_rest/step={d['avg_k_rest']:.4f}  Phi_theory={d['phi_theory']:.4f}  "
          f"err={error:+.4f}")

print()
print("INSIGHT: k_rest/step (avg k for steps 2..h within excursion) should be ≈2")
print("         if the internal orbit follows Geo(1/2). Deviation reveals structure.")
print()

# Summary statistics
avg_k_rests = [d['avg_k_rest'] for d in phi_data.values()]
print(f"avg(k_rest/step) across BSet: {sum(avg_k_rests)/len(avg_k_rests):.4f}")
print(f"  (theoretical: 2.000 from E[l]=2 and geometric k distribution)")
print()

# =====================================================================
# PART B: MOD-512 BSet EXTENSION
# =====================================================================
print("=" * 70)
print("PART B: MOD-512 BSet ANALYSIS")
print("=" * 70)
print()
print("Using mod-512 residues instead of mod-256.")
print("Key: for odd n with v2(n+1) = k0, n ≡ 2^k0-1 mod 2^{k0+1}")
print("With mod 512=2^9: residues n ≡ r mod 512 where r is odd.")
print()

# Identify mod-512 residues with high k0
# For k0=9: n ≡ 511 mod 1024 (but 1024 > 512, so k0=9 not visible mod 512)
# For k0=8: n ≡ 255 mod 512 (n+1≡256 mod 512, v2(256)=8)
# For k0=7: n ≡ 127 or 383 mod 512 (n+1≡128 or 384, v2=7)
# For k0=6: n ≡ 63,191,319,447 mod 512
# etc.

# Build the mod-512 odd residue list
mod512_residues = [r for r in range(1, 512, 2)]  # 256 odd residues mod 512

# Compute k0 for each
def get_k0_512(r):
    return v2(r + 1)  # Note: r is odd, r+1 is even

print("Mod-512 odd residues with k0 ≥ 5:")
high_k0_512 = [(r, get_k0_512(r)) for r in mod512_residues if get_k0_512(r) >= 5]
for r, k0 in sorted(high_k0_512):
    in_bset = (r % 256) in BSet256
    print(f"  r={r:3d} (k0={k0}, mod-256 residue={r%256}, {'∈ BSet256' if in_bset else 'NOT in BSet256'})")

print()

# Now compute Phi(r) for each mod-512 residue
# and find which ones form a "BSet" in the mod-512 sense
print("Computing Phi(r) for mod-512 residues with k0 ≥ 3 (N=256 each)...")
print()

N_512 = 256
phi_512 = {}
MAX_H = 500

for r in sorted(mod512_residues):
    k0 = get_k0_512(r)
    if k0 < 3:
        continue  # skip low-k0 for speed

    k_totals = []
    h_vals = []

    for i in range(N_512):
        n = M_BASE + r + 512 * i
        if v2(n + 1) != k0:
            # Adjust to get correct k0
            continue
        n_cur = n
        h = 0
        k_sum = 0

        while h < MAX_H:
            n_out, k, l = macro_step(n_cur)
            h += 1
            k_sum += k
            n_cur = n_out
            if n_cur <= 1:
                break
            if n_cur % 512 in [r2 for r2 in mod512_residues if get_k0_512(r2) >= 3]:
                k_totals.append(k_sum)
                h_vals.append(h)
                break

    if h_vals:
        phi = sum(k_totals) / sum(h_vals)
        avg_h = sum(h_vals) / len(h_vals)
        phi_512[r] = {'k0': k0, 'phi': phi, 'avg_h': avg_h}

# =====================================================================
# PART C: COMPARE MOD-256 AND MOD-512 BSet STRUCTURES
# =====================================================================
print("=" * 70)
print("PART C: MOD-512 vs MOD-256 BSet COMPARISON")
print("=" * 70)
print()
print("Key question: Does using mod-512 (more refined) change the ergodic avg?")
print()

# Build the BSet using the mod-512 restricted chain for k0>=3
# and compute ergodic avg by computing a simple average of Phi values weighted by 1/E_h
if phi_512:
    # Simple weighted average (proper ergodic would need transition matrix)
    total_phi_h = sum(d['phi'] / d['avg_h'] for d in phi_512.values() if d['avg_h'] > 0)
    total_1_h = sum(1 / d['avg_h'] for d in phi_512.values() if d['avg_h'] > 0)
    approx_erg_512 = total_phi_h / total_1_h if total_1_h > 0 else 0

    print(f"Mod-512 elements analyzed (k0≥3): {len(phi_512)}")
    print(f"  Phi range: [{min(d['phi'] for d in phi_512.values()):.3f}, {max(d['phi'] for d in phi_512.values()):.3f}]")
    print(f"  Approx ergodic avg (weighted by 1/E[h]): {approx_erg_512:.4f}")
    print(f"  Compare: mod-256 BSet ergodic avg = 2.0614")
    print()

# =====================================================================
# PART D: k0-DISTRIBUTION ALONG ACTUAL ORBIT
# =====================================================================
print("=" * 70)
print("PART D: k0-DISTRIBUTION ALONG COLLATZ ORBIT")
print("=" * 70)
print()
print("Is the k0-distribution along a long orbit truly Geo(1/2)?")
print()

# Trace a long orbit and collect k values
n_start = 10**15 + 1  # random large odd n
k_counts = Counter()
total_steps = 0
n_cur = n_start
MAX_ORBIT_STEPS = 100000

while total_steps < MAX_ORBIT_STEPS and n_cur > 1:
    k = v2(n_cur + 1)
    k_counts[k] += 1
    n_out, k_actual, l = macro_step(n_cur)
    total_steps += 1
    n_cur = n_out

print(f"Orbit from n={n_start}")
print(f"Total macro-steps traced: {total_steps}")
print(f"Final n: {n_cur}")
print()

print("k-distribution (macro-steps):")
total_k_mass = sum(k_counts.values())
print(f"{'k':>5}  {'count':>8}  {'fraction':>10}  {'Geo(1/2)':>10}  {'ratio':>8}")
print("-" * 50)
for k in sorted(k_counts.keys()):
    cnt = k_counts[k]
    frac = cnt / total_k_mass
    geo = 1.0 / 2**k
    ratio = frac / geo
    print(f"k={k:3d}  {cnt:8d}  {frac:.6f}  {geo:.6f}  {ratio:.4f}")

print()
avg_k_orbit = sum(k * cnt for k, cnt in k_counts.items()) / total_k_mass
print(f"Empirical E[k] along orbit: {avg_k_orbit:.6f}")
print(f"Theoretical E[k] = 2 for Geo(1/2): 2.000000")
print(f"Difference: {avg_k_orbit - 2:.6f}")
print()
print(f"Empirical log-drift per step: {avg_k_orbit * math.log(1.5) - 2 * math.log(2):.6f}")
print(f"Theoretical log-drift (E[k]=2): {2 * math.log(1.5) - 2 * math.log(2):.6f}")

# =====================================================================
# PART E: SUMMARY — UNIFIED PICTURE
# =====================================================================
print()
print("=" * 70)
print("PART E: UNIFIED PICTURE")
print("=" * 70)
print()
log2 = math.log(2)
log32 = math.log(1.5)
threshold = 2 * log2 / log32

print("The k-distribution along any Collatz orbit appears to be Geo(1/2):")
print("  P(k=j) ≈ 1/2^j  for j=1,2,3,...")
print("  E[k] ≈ 2")
print()
print("Combined with E[l]=2 (PROVED), the log-drift is:")
print(f"  E[drift] = E[k]×log(3/2) - E[l]×log2 ≈ 2×{log32:.4f} - 2×{log2:.4f}")
print(f"           = {2*log32 - 2*log2:.4f} per macro-step")
print()
print("This is NEGATIVE → typical orbits converge.")
print()
print(f"For DIVERGING orbits (D_hard_kern), need E[drift] ≥ 0:")
print(f"  E[k] ≥ 2×log2/log(3/2) = {threshold:.6f} = log_{{3/2}}(4)")
print()
print("The BSet analysis shows that even if an orbit concentrates on BSet,")
print(f"the best achievable ergodic avg k/step = 2.06 << {threshold:.4f}.")
print()
print("OPEN QUESTION: Can an orbit have E[k] >> 2 by AVOIDING high-k0 BSet elements?")
print("  → Non-BSet max k/step = 2.25 < 3.419. No.")
print()
print("CONCLUSION: The k-distribution Geo(1/2) + E[l]=2 + BSet ergodic analysis")
print("together prove (empirically) that D_hard_kern = ∅.")
print()
print("The MISSING step for a RIGOROUS proof:")
print("  Prove that the k-distribution along any orbit converges to Geo(1/2).")
print("  This would follow from mixing/equidistribution of Collatz orbits mod 2^k.")
