"""
92_bset_threshold_all_residues.py
==================================
For EACH of the 128 odd residues mod 256:
  - Compute long-run avg k/step using N large-n starting points
  - Classify as BSet or non-BSet
  - Find the EXACT threshold in avg k/step that separates the two sets

Also: show the exact arithmetic formula for E[l] (universal constant = 2.0)
and compute it for ALL 128 residues (not just the 15 BSet elements).

Key question: Is "avg k/step > threshold" the EXACT criterion for BSet membership?
If so, what is the threshold?
"""
import sys, time
import numpy as np
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
BList = sorted(BSet)
ODD_RES = [r for r in range(1, 256, 2)]  # 128 odd residues

M_BASE = 10**12
M_PERIOD = 256
MAX_H = 200

print("=" * 70)
print("PART 1: E[l] FOR ALL 128 ODD RESIDUES (EXACT 256-POINT COMPUTATION)")
print("=" * 70)
print()
print("l = v₂(3^k₀ × m - 1) for k₀ = v₂(r+1) and 256 odd m-values.")
print()

t0 = time.time()
E_l_by_res = {}
k0_by_res = {}
for r in ODD_RES:
    k0 = v2(r+1)
    k0_by_res[r] = k0
    pow3k0 = 3**k0
    l_vals = []
    for m_idx in range(M_PERIOD):
        m = 2*m_idx + 1
        val = pow3k0 * m - 1
        l_vals.append(v2(val))
    E_l_by_res[r] = sum(l_vals) / len(l_vals)

# Check: is E[l] universally ≈ 2?
E_l_vals = list(E_l_by_res.values())
print(f"E[l] statistics over all 128 odd residues:")
print(f"  min = {min(E_l_vals):.4f}")
print(f"  max = {max(E_l_vals):.4f}")
print(f"  mean = {sum(E_l_vals)/len(E_l_vals):.4f}")
print(f"  std = {np.std(E_l_vals):.6f}")

# Show any outliers (E[l] != 2.0 ± 0.01)
outliers = [(r, E_l_by_res[r]) for r in ODD_RES if abs(E_l_by_res[r]-2.0) > 0.02]
if outliers:
    print(f"\nResidues with E[l] far from 2.0:")
    for r, el in sorted(outliers, key=lambda x: abs(x[1]-2)):
        k0 = k0_by_res[r]
        print(f"  r={r:3d} k0={k0}: E[l]={el:.4f}")
else:
    print(f"\nAll 128 residues have E[l] within ±0.02 of 2.0. E[l]=2 is UNIVERSAL.")

# =====================================================================
# PART 2: LONG-RUN AVG k/STEP FOR ALL 128 RESIDUES
# =====================================================================
print()
print("=" * 70)
print("PART 2: LONG-RUN AVG k/STEP AND BSet THRESHOLD")
print("=" * 70)
print()
print(f"Computing long-run avg k/step for all 128 odd residues...")
print(f"Using {M_PERIOD} large-n starting points per residue (m ~ 10^12)")
print()

t1 = time.time()
avg_k_step = {}  # r -> avg k per macro-step along the orbit to BSet hit
avg_h = {}       # r -> avg hop count to BSet hit
convergence_rate = {}  # r -> fraction that converge before BSet hit

for r in ODD_RES:
    k0 = k0_by_res[r]
    step = 2**(k0 + 1)

    k_sums = []
    h_vals = []
    n_converged = 0

    for m_idx in range(M_PERIOD):
        m = M_BASE + 2*m_idx + 1
        n = step * m - 1

        # Fix if k != k0
        if v2(n+1) != k0:
            n = ((M_BASE // step) + 1) * step + r
            if n < M_BASE:
                n += step * ((M_BASE - n) // step + 1)

        n_cur = n
        h = 0
        k_sum = 0
        outcome = None

        while h < MAX_H:
            n_cur, k, l = macro_step(n_cur)
            h += 1
            k_sum += k
            if n_cur <= 1:
                n_converged += 1
                outcome = 'converged'
                break
            if n_cur % 256 in BSet:
                outcome = 'bset'
                break

        if outcome == 'bset':
            k_sums.append(k_sum)
            h_vals.append(h)

    if h_vals:
        avg_k_step[r] = sum(k_sums) / sum(h_vals)
        avg_h[r] = sum(h_vals) / len(h_vals)
        convergence_rate[r] = n_converged / M_PERIOD
    else:
        avg_k_step[r] = None  # all converged
        avg_h[r] = None
        convergence_rate[r] = 1.0

print(f"[Done in {time.time()-t1:.1f}s]\n")

# Sort by avg k/step
sorted_res = sorted([r for r in ODD_RES if avg_k_step[r] is not None],
                    key=lambda r: avg_k_step[r])

# Show the BSet/non-BSet boundary
print("=== AVG k/STEP SORTED — BSet/non-BSet BOUNDARY ===\n")
print(f"{'rank':>5}  {'r':>5}  {'k0':>4}  {'avg_k/step':>12}  {'avg_h':>7}  {'conv%':>7}  {'BSet?':>8}")
print("-" * 65)

bset_k_steps = [avg_k_step[r] for r in ODD_RES if r in BSet and avg_k_step[r] is not None]
non_bset_k_steps = [avg_k_step[r] for r in ODD_RES if r not in BSet and avg_k_step[r] is not None]

if bset_k_steps and non_bset_k_steps:
    min_bset = min(bset_k_steps)
    max_non_bset = max(non_bset_k_steps)
    threshold = (min_bset + max_non_bset) / 2

    # Show all residues within 0.2 of the threshold
    for i, r in enumerate(sorted_res):
        k0 = k0_by_res[r]
        ks = avg_k_step[r]
        h = avg_h[r]
        cr = convergence_rate[r]
        in_bset = r in BSet
        if abs(ks - threshold) < 0.2 or (i < 5) or (i >= len(sorted_res)-5):
            marker = " ← BSet boundary" if abs(ks - threshold) < 0.05 else ""
            print(f"  {i+1:3d}  r={r:3d}  k={k0}  avg_k/step={ks:.4f}  h={h:.2f}  conv={100*cr:.1f}%  {'BSet ✓' if in_bset else 'non-BSet':>8}{marker}")

print()
print(f"BSet avg k/step range:     [{min(bset_k_steps):.4f}, {max(bset_k_steps):.4f}]")
print(f"Non-BSet avg k/step range: [{min(non_bset_k_steps):.4f}, {max(non_bset_k_steps):.4f}]")
print(f"Gap (min_BSet - max_nonBSet): {min(bset_k_steps) - max(non_bset_k_steps):.4f}")
print(f"Estimated threshold: {threshold:.4f}")
print()
if min(bset_k_steps) > max(non_bset_k_steps):
    print("✓ CLEAN SEPARATION: BSet and non-BSet are PERFECTLY SEPARATED by avg k/step!")
else:
    print("✗ OVERLAP: BSet and non-BSet overlap in avg k/step range.")

# =====================================================================
# PART 3: EXACT BSET CLASSIFICATION SORTED BY k/step
# =====================================================================
print()
print("=== FULL CLASSIFICATION: 128 ODD RESIDUES BY AVG k/STEP ===\n")
print(f"{'r':>5}  {'k0':>4}  {'avg_k/step':>12}  {'avg_h':>8}  {'conv%':>7}  {'BSet?':>8}")
print("-" * 60)

# Show only the boundary region: 20 below and 20 above threshold
boundary_res = [(r, avg_k_step[r]) for r in sorted_res if avg_k_step[r] is not None]
# Find threshold index
idx_thresh = next(i for i, (r, ks) in enumerate(boundary_res) if ks >= threshold)

print("... (non-BSet residues with low avg k/step) ...")
for i in range(max(0, idx_thresh-10), len(boundary_res)):
    r, ks = boundary_res[i]
    k0 = k0_by_res[r]
    h = avg_h[r]
    cr = convergence_rate[r]
    in_bset = r in BSet
    arrow = " ← THRESHOLD" if idx_thresh-1 <= i <= idx_thresh else ""
    print(f"  r={r:3d}  k={k0}  avg_k/step={ks:.4f}  h={h:.2f}  conv={100*cr:.1f}%  {'BSet ✓' if in_bset else 'non-BSet':>8}{arrow}")
    if i > idx_thresh + 10: break

# =====================================================================
# PART 4: THE SINGLE-STEP DRIFT FORMULA
# =====================================================================
print()
print("=" * 70)
print("PART 4: SINGLE-STEP DRIFT = k₀ × log₂3 - E[l] FOR ALL RESIDUES")
print("=" * 70)
print()
import math
log23 = math.log2(3)  # = 1.58496...

print("Single-step log₂-drift = k₀ × log₂3 - E[l]  (using E[l]≈2.0)")
print()

drift_by_k0 = defaultdict(list)
for r in ODD_RES:
    k0 = k0_by_res[r]
    drift = k0 * log23 - E_l_by_res[r]
    drift_by_k0[k0].append((r, drift, r in BSet))

print(f"{'k0':>4}  {'single-step drift':>18}  {'n_bset':>8}  {'n_total':>8}")
for k0 in sorted(drift_by_k0.keys()):
    items = drift_by_k0[k0]
    drifts = [d for _,d,_ in items]
    n_bset = sum(1 for _,_,b in items if b)
    n_total = len(items)
    print(f"  k0={k0}: drift={drifts[0]:+.4f}  n_bset={n_bset}/{n_total}")

print()
print(f"Drift crosses zero at k0×log₂3 > 2  →  k0 > 2/log₂3 = 2/1.585 = {2/log23:.3f}")
print(f"→ k0 ≥ 2 gives positive single-step drift")
print(f"→ But BSet requires avg k/step > threshold ≈ {threshold:.3f}")
print()
print(f"OBSERVATION: BSet threshold (avg k/step ≈ {threshold:.3f}) is NOT the same as")
print(f"single-step drift threshold (k0 ≥ 2). The BSet condition captures MULTI-STEP behavior.")

# =====================================================================
# PART 5: WHY BSet THRESHOLD ≈ 2.15-2.20?
# =====================================================================
print()
print("=" * 70)
print("PART 5: BSet THRESHOLD AND D_hard_kern GAP")
print("=" * 70)
print()
print(f"BSet threshold (from data): ≈ {threshold:.4f}")
print(f"D_hard_kern threshold (from Theorem 179): 3.419")
print(f"Max observed avg k/step (r=255): 2.279")
print()
print(f"THREE THRESHOLDS:")
print(f"  1. BSet membership:    avg k/step ≥ {threshold:.3f}")
print(f"  2. Orbit convergence:  avg k/step <  {threshold:.3f}  (converges to 1)")
print(f"  3. D_hard_kern:        avg k/step ≥  3.419           (diverges)")
print()
print(f"'Safe zone' for BSet elements: avg k/step in [{threshold:.3f}, 2.279]")
print(f"D_hard_kern threshold (3.419) is {3.419-2.279:.3f} above the maximum BSet avg k/step.")
print()
print(f"ALL 15 BSet elements have avg k/step in [{min(bset_k_steps):.3f}, {max(bset_k_steps):.3f}] — WELL below 3.419.")
print()

# Show BSet elements sorted by avg k/step
print("BSet elements sorted by avg k/step:")
bset_sorted = [(r, avg_k_step[r], k0_by_res[r]) for r in BList if avg_k_step[r] is not None]
bset_sorted.sort(key=lambda x: x[1])
for r, ks, k0 in bset_sorted:
    h = avg_h[r]
    print(f"  r={r:3d} (k0={k0}): avg_k/step={ks:.4f}  avg_h={h:.2f}")
print()
print(f"Max BSet avg k/step: {max(bset_k_steps):.4f}")
print(f"D_hard_kern gap from max BSet: {3.419 - max(bset_k_steps):.4f}")

# Compare with non-BSet near-boundary
print()
print("Non-BSet near-boundary elements (closest to BSet):")
near_bset = [(r, avg_k_step[r], k0_by_res[r]) for r in ODD_RES
             if r not in BSet and avg_k_step[r] is not None and avg_k_step[r] > threshold - 0.15]
near_bset.sort(key=lambda x: -x[1])
for r, ks, k0 in near_bset[:10]:
    h = avg_h[r]
    cr = convergence_rate[r]
    print(f"  r={r:3d} (k0={k0}): avg_k/step={ks:.4f}  avg_h={h:.2f}  conv%={100*cr:.1f}%")
