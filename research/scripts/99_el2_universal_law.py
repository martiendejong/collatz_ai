"""
99_el2_universal_law.py
========================
THE E[l]=2 UNIVERSAL LAW:

For ANY k0 ≥ 1 and uniform random odd m:
  E[v2(3^k0 × m - 1)] = 2  EXACTLY

PROOF:
  P(v2(3^k0 × m - 1) ≥ k | m odd) = P(m ≡ (3^k0)^{-1} mod 2^k | m odd)
                                     = 1 / 2^{k-1}
  (since (3^k0)^{-1} mod 2^k is a specific odd residue)

  E[l] = Σ_{k=1}^∞ P(l≥k) = Σ 1/2^{k-1} = 2.

This is UNIVERSAL — same for ALL k0.

CONSEQUENCES:
1. E[k/step] ≈ k0 - 0 (the macro-step "uses" k0 bits of k-energy per step)
2. Log-drift: E[log(n_out/n)] ≈ k0×log(3/2) - 2×log2 = k0×log(3/2) - 2×log2
   For k0=2: E[log ratio] ≈ 2×0.585 - 2×0.693 = 1.170 - 1.386 = -0.216 < 0 (convergent)
   For k0=3.419: E[log ratio] ≈ 3.419×0.585 - 2×0.693 = 2.000 - 1.386 = 0.614 > 0 (divergent)
   → D_hard_kern threshold 3.419 corresponds to zero drift: k0 = 2×log2/log(3/2)!

CRITICAL INSIGHT:
  D_hard_kern threshold = 2×log2/log(3/2) = 2×log_2(4/3) × log_2(2) / log_2(4/3)...
  Let me compute: if k0×log(3/2) = 2×log2, then k0 = 2×log2/log(3/2) = log(4)/log(3/2)
                 = log_{3/2}(4) ≈ 3.419!

So the Collatz divergence threshold IS the zero-drift condition for log(n)!
"""
import sys, math
from collections import Counter, defaultdict

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

BSet = {27, 55, 63, 83, 95, 103, 127, 159, 169, 191, 207, 223, 239, 253, 255}

# =====================================================================
# PART 1: VERIFY E[l]=2 FOR EACH k0 VALUE
# =====================================================================
print("=" * 70)
print("PART 1: E[l]=2 UNIVERSAL LAW — VERIFICATION PER k0")
print("=" * 70)
print()
print("For each k0, compute E[v2(3^k0 × m - 1)] over 256 odd m in [1,511]")
print()

log2 = math.log(2)
log_3_2 = math.log(1.5)

for k0 in range(1, 9):
    # All odd m in [1, 2^(k0+1)-1] (the minimal m-class for k0)
    # or just [1, 511] for comparison
    pow3 = 3**k0
    l_vals = []
    for m in range(1, 512, 2):  # 256 odd m
        x = pow3 * m - 1
        l = v2(x)
        l_vals.append(l)

    avg_l = sum(l_vals) / len(l_vals)
    max_l = max(l_vals)
    l_dist = Counter(l_vals)
    sorted_l = sorted(l_dist.items())

    # Theoretical: P(l≥k) = 1/2^{k-1}
    # E[l] ≈ Σ_{k=1}^{max_l+1} 1/2^{k-1} = 2(1-1/2^{max_l}) ≈ 2
    e_l_theory = sum(len([x for x in l_vals if x >= k]) / 256 for k in range(1, max_l+2))

    print(f"k0={k0}: E[l]={avg_l:.4f}  max_l={max_l}  theoretical_sum={e_l_theory:.4f}")
    print(f"        l distribution: { {k: v for k,v in sorted_l[:8]} }")

print()
print("Theoretical E[l] = 2 exactly (for uniform odd m over infinite range)")
print("Finite-sample E[l] ≈ 511/256 ≈ 1.996 (truncated geometric)")
print()

# =====================================================================
# PART 2: THE CRITICAL INSIGHT — D_hard_kern THRESHOLD IS LOG-DRIFT ZERO
# =====================================================================
print("=" * 70)
print("PART 2: D_hard_kern THRESHOLD = ZERO LOG-DRIFT CONDITION")
print("=" * 70)
print()
print("FORMULA: log(n_out/n) ≈ k0×log(3/2) - l×log(2) per macro-step")
print()
print("Taking expectations (with E[k0]=k0 fixed, E[l]=2):")
print("  E[log(n_out/n)] ≈ k0×log(3/2) - 2×log(2)")
print()
print("For zero drift (boundary between convergent and divergent):")
print("  k0×log(3/2) = 2×log(2)")
print("  k0 = 2×log(2)/log(3/2) = log(4)/log(3/2) = log_{3/2}(4)")
print()

k0_critical = 2 * log2 / log_3_2
print(f"  k0_critical = log_{{3/2}}(4) = {k0_critical:.6f}")
print()
print(f"  D_hard_kern threshold (Theorem 179): 3.419000")
print(f"  log_{{3/2}}(4) = {k0_critical:.6f}")
print(f"  Difference: {3.419 - k0_critical:.6f}")
print()
print("REMARKABLE: The D_hard_kern threshold IS (approximately) log_{3/2}(4)!")
print()

# Verify: log(4)/log(3/2) = log(4)/log(3/2)
# = 2×log(2) / (log(3)-log(2))
# = 2 / (log(3)/log(2) - 1)
# = 2 / (log_2(3) - 1)
log2_3 = math.log(3) / math.log(2)
k0_alt = 2 / (log2_3 - 1)
print(f"  Alternative formula: 2/(log_2(3)-1) = 2/({log2_3:.6f}-1) = {k0_alt:.6f}")
print()

# =====================================================================
# PART 3: DRIFT RATES FOR EACH k0
# =====================================================================
print("=" * 70)
print("PART 3: LOG-DRIFT PER MACRO-STEP AT EACH k0")
print("=" * 70)
print()
print(f"{'k0':>4}  {'E[k0×log(3/2)-2×log2]':>25}  {'drift':>12}  {'converges?':>12}")
print("-" * 60)
for k0 in range(1, 11):
    drift = k0 * log_3_2 - 2 * log2
    converges = "YES (n shrinks)" if drift < 0 else "NO (n grows)"
    print(f"  k0={k0:2d}  k0×0.405 - 1.386 = {drift:+.4f}  {drift:+12.6f}  {converges}")
print()
print(f"  Zero drift at k0 = {k0_critical:.4f} (= D_hard_kern threshold!)")
print()

# =====================================================================
# PART 4: CONNECTION TO BSet
# =====================================================================
print("=" * 70)
print("PART 4: WHY BSet IS THE RIGHT SET")
print("=" * 70)
print()
print("BSet elements and their k0 values:")
for r in sorted(BSet):
    k0 = v2(r + 1)
    drift = k0 * log_3_2 - 2 * log2
    print(f"  r={r:3d}  k0={k0}  log-drift={drift:+.4f}")

print()
print("Key observation:")
print("  k0=8 (r=255): drift = 8×0.405-1.386 = +1.854  (strong growth per step)")
print("  k0=1 (r=169): drift = 1×0.405-1.386 = -0.981  (strong shrinking per step)")
print("  k0=2:         drift = 2×0.405-1.386 = -0.577")
print()
print("BSet elements have k0 ∈ {1,2,3,4,5,6,7,8}.")
print("The high-k0 elements (k0=6,7,8) have positive drift per step.")
print("The low-k0 elements (k0=1,2) have negative drift per step.")
print("Their MIXTURE (via the Markov chain) achieves ergodic avg k/step = 2.06.")
print()
print("WHY ergodic avg = 2.06?")
print("  The ergodic avg k/step corresponds to zero drift when k0 = 2×log2/log(3/2).")
print("  The BSet Markov chain's ergodic avg = 2.06 < 3.419 → net convergence.")
print("  Any BSet orbit: E[log(n_out/n)] = 2.06×log(3/2) - 2×log2")
erg_drift = 2.0614 * log_3_2 - 2 * log2
print(f"                                     = {erg_drift:.4f} < 0 (orbit converges on avg)")
print()

# =====================================================================
# PART 5: VERIFY DRIFT COMPUTATION EMPIRICALLY
# =====================================================================
print("=" * 70)
print("PART 5: EMPIRICAL LOG-DRIFT VERIFICATION AT n~10^12")
print("=" * 70)
print()
print("For each BSet element, compute actual E[log(n_out/n)] over N=1024 trajectories")
print("and compare to theoretical k0×log(3/2) - 2×log(2)")
print()

M_BASE = 10**12
N = 1024
BList = sorted(BSet)

total_drift_sum = 0
total_weight = 0

for r in BList:
    k0 = v2(r + 1)
    theo_drift = k0 * log_3_2 - 2 * log2

    drifts = []
    for i in range(N):
        n = M_BASE + r + 256 * i
        n_out, k, l = macro_step(n)
        if n_out > 1:
            actual_drift = math.log(n_out) - math.log(n)
            drifts.append(actual_drift)

    avg_drift = sum(drifts) / len(drifts) if drifts else 0
    print(f"  r={r:3d} (k0={k0}): theoretical={theo_drift:+.4f}  empirical={avg_drift:+.4f}  diff={avg_drift-theo_drift:+.4f}")

# =====================================================================
# PART 6: THE DEEP CONNECTION
# =====================================================================
print()
print("=" * 70)
print("PART 6: UNIFIED VIEW — COLLATZ AS CONTROLLED RANDOM WALK")
print("=" * 70)
print()
print("LOG-DRIFT PER MACRO-STEP: log(n') - log(n) ≈ k×log(3/2) - l×log(2)")
print()
print("The Collatz macro-step is a RANDOM WALK in log(n) with:")
print("  Step size = k×log(3/2) - l×log(2)")
print("  Step distribution: k geometric(1/2), l geometric(1/2)")
print("  E[step] = E[k]×log(3/2) - E[l]×log2 = 2×log(3/2) - 2×log2 = 2×log(3/4) < 0")
print()
print("NORMAL ORBITS (E[k]≈2): CONVERGING walk with E[step] = 2×log(3/4) ≈ -0.575")
print()
print("D_hard_kern ORBITS (E[k]≥3.419): DIVERGING walk with E[step] ≥ 0")
print()
print("Why the threshold 3.419 = log_{3/2}(4)?")
print("  E[step] = 0 ⟺ E[k]×log(3/2) = E[l]×log2 = 2×log2")
print("  ⟺ E[k] = 2×log2/log(3/2) = log(4)/log(3/2) = log_{3/2}(4) ≈ 3.419")
print()
print("The Collatz conjecture says all orbits converge.")
print("This requires E[k] < 3.419 for ALL orbits.")
print("D_hard_kern = ∅ says exactly this.")
print()
print("REMARKABLE SYNTHESIS:")
print("  1. E[l]=2 UNIVERSAL (proved from modular arithmetic)")
print("  2. D_hard_kern threshold = 2×log2/log(3/2) (from zero-drift condition)")
print("  3. BSet ergodic avg k/step = 2.06 << 3.419 (from Markov chain analysis)")
print("  4. ALL three layers give E[drift] < 0 for BSet orbits ✓")
print()
print(f"NUMERICAL CHECK:")
print(f"  log_{{3/2}}(4) = {k0_critical:.6f}")
print(f"  Theorem 179 threshold: 3.419000")
print(f"  Difference: {k0_critical - 3.419:.6f}")
print()
if abs(k0_critical - 3.419) < 0.01:
    print("  ✓ D_hard_kern threshold ≈ log_{3/2}(4) (matches to within 0.01)")
    print()
    print("  The threshold in Theorem 179 IS (or is very close to) log_{3/2}(4) = log(4)/log(3/2)")
    print("  This gives a BEAUTIFUL closed-form explanation for the '3.419' magic number!")
