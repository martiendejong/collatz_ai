"""
109_full_mod256_markov.py
==========================
FULL MOD-256 MARKOV CHAIN: ALL 128 ODD RESIDUES

For each odd residue r in [1,255], compute the transition probabilities
P(r -> r') under the MACRO-STEP map (not restricted to BSet).

Key questions:
1. Is the chain ergodic (irreducible + aperiodic)?
2. What is the stationary distribution? Is it uniform?
3. What is the spectral gap?
4. Does the stationary distribution weight BSet elements more or less?

Connection to equidistribution: if Collatz equidistribution holds mod 256,
then the stationary distribution of this chain IS the uniform distribution.
How far is it from uniform?
"""
import sys, math
from collections import Counter

sys.stdout.reconfigure(encoding='utf-8')

def v2(x):
    if x == 0: return 999
    c = 0
    while x % 2 == 0: x //= 2; c += 1
    return c

def macro_step(n):
    k = v2(n + 1); m = (n + 1) >> k; x = m * (3**k) - 1; l = v2(x); return x >> l, k, l

BSet = {27,55,63,83,95,103,127,159,169,191,207,223,239,253,255}
odd_res = list(range(1, 256, 2))  # 128 odd residues
N = len(odd_res)  # = 128
idx = {r: i for i, r in enumerate(odd_res)}

# =====================================================================
# PART 1: COMPUTE TRANSITION MATRIX P[r][r'] MOD 256
# For each residue r, sample n = r + 256k (k=0,...,N_SAMP-1)
# and track where n' lands mod 256.
# =====================================================================
print("=" * 70)
print("PART 1: TRANSITION MATRIX (128x128) ON ODD RESIDUES MOD 256")
print("=" * 70)
print()

N_SAMP = 512  # per residue
P = [[0.0] * N for _ in range(N)]

for i, r in enumerate(odd_res):
    counts = Counter()
    for k in range(N_SAMP):
        n = r + 256 * k
        n_out, _, _ = macro_step(n)
        r_out = n_out % 256
        if r_out % 2 == 1 and 1 <= r_out <= 255:
            counts[r_out] += 1
    total = sum(counts.values())
    for r2, c in counts.items():
        if r2 in idx:
            P[i][idx[r2]] = c / total

print(f"Matrix computed: {N}x{N}, each row from {N_SAMP} samples.")
print()

# =====================================================================
# PART 2: STATIONARY DISTRIBUTION (POWER ITERATION)
# =====================================================================
print("=" * 70)
print("PART 2: STATIONARY DISTRIBUTION")
print("=" * 70)
print()

pi = [1.0/N] * N

for iteration in range(2000):
    pi_new = [0.0] * N
    for i in range(N):
        for j in range(N):
            pi_new[j] += pi[i] * P[i][j]
    max_diff = max(abs(pi_new[j] - pi[j]) for j in range(N))
    pi = pi_new
    if max_diff < 1e-10:
        print(f"Converged at iteration {iteration+1}")
        break

# Compare to uniform
uniform = 1.0/N
print()
print(f"Uniform = 1/128 = {uniform:.6f}")
print()

# L1 deviation from uniform
l1_dev = sum(abs(pi[i] - uniform) for i in range(N))
max_dev = max(abs(pi[i] - uniform) for i in range(N))
print(f"L1 deviation from uniform: {l1_dev:.6f}")
print(f"Max deviation from uniform: {max_dev:.6f} ({max_dev/uniform*100:.1f}%)")
print()

# Most visited residues
sorted_pi = sorted(enumerate(pi), key=lambda x: -x[1])
print("Top 15 most visited odd residues:")
for rank, (i, p) in enumerate(sorted_pi[:15]):
    r = odd_res[i]
    is_bset = "BSet*" if r in BSet else ""
    print(f"  #{rank+1:2d}: r={r:3d} (k0={v2(r+1)}), pi={p:.6f} ({p/uniform*100:.1f}% of uniform) {is_bset}")

print()
print("Bottom 5 least visited:")
for rank, (i, p) in enumerate(sorted_pi[-5:]):
    r = odd_res[i]
    is_bset = "BSet*" if r in BSet else ""
    print(f"  #{N-5+rank+1:3d}: r={r:3d} (k0={v2(r+1)}), pi={p:.6f} ({p/uniform*100:.1f}%) {is_bset}")

# =====================================================================
# PART 3: STATIONARY WEIGHT OF BSet vs non-BSet
# =====================================================================
print()
print("=" * 70)
print("PART 3: BSet vs NON-BSet STATIONARY WEIGHT")
print("=" * 70)
print()

bset_weight = sum(pi[idx[r]] for r in BSet if r in idx)
nonbset_weight = 1.0 - bset_weight

print(f"BSet elements (15): stationary weight = {bset_weight:.6f}")
print(f"Non-BSet elements (113): stationary weight = {nonbset_weight:.6f}")
print()
print(f"Uniform BSet weight = 15/128 = {15/128:.6f}")
print(f"BSet overweight factor = {bset_weight/(15/128):.4f}x")
print()

# Per BSet element
print("Stationary weights for BSet elements:")
for r in sorted(BSet):
    K = v2(r+1)
    p = pi[idx[r]]
    print(f"  r={r:3d} (K={K}): pi={p:.6f} ({p/(1/128)*100:.1f}% of uniform)")

# =====================================================================
# PART 4: ERGODIC AVERAGE k PER STEP (OVERALL ORBIT)
# =====================================================================
print()
print("=" * 70)
print("PART 4: ERGODIC AVERAGE k PER MACRO-STEP (FULL ORBIT)")
print("=" * 70)
print()

# Compute avg k for each residue (the k0 of the macro-step starting from r)
# = E[v2(n+1)] where n ≡ r mod 256
avg_k_from_r = {}
for r in odd_res:
    k_vals = []
    for k_iter in range(N_SAMP):
        n = r + 256 * k_iter
        k_vals.append(v2(n+1))
    avg_k_from_r[r] = sum(k_vals)/len(k_vals)

ergodic_k = sum(pi[idx[r]] * avg_k_from_r[r] for r in odd_res)
print(f"Ergodic avg k0 per step (from stationary pi) = {ergodic_k:.6f}")
print(f"Threshold                                     = {2*math.log(2)/math.log(1.5):.6f}")
print(f"Gap                                           = {2*math.log(2)/math.log(1.5)-ergodic_k:.6f}")
print()
print("By k0 class:")
for k_val in range(1, 10):
    k_residues = [r for r in odd_res if avg_k_from_r[r] == k_val]
    if k_residues:
        total_weight = sum(pi[idx[r]] for r in k_residues)
        n_res = len(k_residues)
        print(f"  k0={k_val}: {n_res:3d} residues, total pi = {total_weight:.5f}, "
              f"contribution = {total_weight*k_val:.5f}")

# =====================================================================
# PART 5: SPECTRAL GAP ANALYSIS
# =====================================================================
print()
print("=" * 70)
print("PART 5: SPECTRAL ANALYSIS AND MIXING TIME")
print("=" * 70)
print()

try:
    import numpy as np
    M_np = np.array(P)
    eigenvalues = np.linalg.eigvals(M_np)
    eigenvalues_real = sorted(eigenvalues.real, reverse=True)
    lambda2 = eigenvalues_real[1]
    spectral_gap = 1.0 - lambda2
    print(f"Leading eigenvalue: {eigenvalues_real[0]:.6f} (should be 1.000)")
    print(f"Second eigenvalue: {lambda2:.6f}")
    print(f"Spectral gap = 1 - lambda_2 = {spectral_gap:.6f}")
    print(f"Mixing time estimate ~ 1/gap = {1/spectral_gap:.2f} macro-steps")
    print()
    print("Top 5 eigenvalues:")
    for i, ev in enumerate(eigenvalues_real[:5]):
        print(f"  lambda_{i+1} = {ev:.6f}")
except ImportError:
    print("NumPy not available; skipping spectral analysis.")
except Exception as e:
    print(f"Spectral analysis error: {e}")

# =====================================================================
# PART 6: CHECK IF STATIONARY = UNIFORM (DETAILED)
# =====================================================================
print()
print("=" * 70)
print("PART 6: TEST UNIFORMITY HYPOTHESIS")
print("=" * 70)
print()
print("HYPOTHESIS: Collatz equidistribution => stationary pi = uniform")
print()
print("L1 deviation test:")
print(f"  Observed L1 deviation = {l1_dev:.6f}")
print(f"  For truly uniform with N={N} samples of {N_SAMP}:")
print(f"    Expected sampling noise ~ sqrt(N/N_SAMP) = {(N/N_SAMP)**0.5:.4f}")
print()
print("Residues MOST above uniform:")
for i, p in sorted_pi[:8]:
    r = odd_res[i]
    excess = (p - uniform)/uniform*100
    print(f"  r={r:3d}: pi={p:.6f} (+{excess:.1f}% above uniform), k0={v2(r+1)}, BSet={r in BSet}")
print()
print("Residues MOST below uniform:")
for i, p in sorted_pi[-8:]:
    r = odd_res[i]
    deficit = (p - uniform)/uniform*100
    print(f"  r={r:3d}: pi={p:.6f} ({deficit:.1f}% below uniform), k0={v2(r+1)}, BSet={r in BSet}")
