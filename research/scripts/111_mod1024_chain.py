"""
111_mod1024_chain.py
=====================
MOD-1024 MARKOV CHAIN — 512 ODD RESIDUES
Third data point for the EXPANDER CONJECTURE:
  spectral_gap(Collatz mod 2^N) >= 0.90 for all N >= 8.

Data so far:
  Mod-256 (N=128 states):  gap = 0.938189
  Mod-512 (N=256 states):  gap = 0.920260
  Mod-1024 (N=512 states): gap = ???

Also investigate: does the second eigenvalue pattern continue?
"""
import sys, math, numpy as np
from collections import Counter

sys.stdout.reconfigure(encoding='utf-8')

def v2(x):
    if x == 0: return 999
    c = 0
    while x % 2 == 0: x //= 2; c += 1
    return c

def macro_step(n):
    k = v2(n + 1); m = (n + 1) >> k; x = m * (3**k) - 1; l = v2(x)
    return x >> l, k, l

BSet_256 = {27,55,63,83,95,103,127,159,169,191,207,223,239,253,255}

odd_res_1024 = list(range(1, 1024, 2))  # 512 odd residues
N1024 = 512
idx1024 = {r: i for i, r in enumerate(odd_res_1024)}

# BSet mod 1024: each mod-256 BSet element gives 4 mod-1024 residues
BSet_1024 = set()
for r256 in BSet_256:
    for k in range(4):
        r1024 = r256 + 256 * k
        if r1024 % 2 == 1 and r1024 < 1024:
            BSet_1024.add(r1024)

print(f"States: {N1024} (odd residues mod 1024)")
print(f"BSet_1024: {len(BSet_1024)} elements (expected 60 = 15 × 4)")
print()

N_SAMP = 256  # per residue — enough for spectral gap estimate

print(f"Computing {N1024}x{N1024} transition matrix ({N_SAMP} samples/state)...")
P1024 = np.zeros((N1024, N1024))
for i, r in enumerate(odd_res_1024):
    if i % 64 == 0:
        print(f"  Progress: {i}/{N1024} states done...")
    counts = np.zeros(N1024)
    for k in range(N_SAMP):
        n = r + 1024 * k
        n_out, _, _ = macro_step(n)
        r_out = n_out % 1024
        if r_out % 2 == 1 and r_out in idx1024:
            counts[idx1024[r_out]] += 1
    total = counts.sum()
    if total > 0:
        P1024[i] = counts / total

print("Matrix computed.")
print()

# Spectral analysis
print("Computing eigenvalues (this may take a moment for 512x512)...")
vals1024 = np.linalg.eigvals(P1024)
vals1024_sorted = sorted(vals1024.real, reverse=True)
gap1024 = 1 - vals1024_sorted[1]

print("Top 10 eigenvalues (mod-1024):")
for i in range(10):
    print(f"  lambda_{i+1} = {vals1024_sorted[i]:.8f}")
print()
print(f"Spectral gap (mod-1024) = {gap1024:.8f}")
print()

# Stationary distribution (power iteration)
pi1024 = np.ones(N1024) / N1024
for _ in range(200):
    pi1024_new = pi1024 @ P1024
    if np.max(np.abs(pi1024_new - pi1024)) < 1e-8:
        break
    pi1024 = pi1024_new

uniform1024 = 1.0 / N1024
l1_1024 = np.sum(np.abs(pi1024 - uniform1024))
max_dev_1024 = np.max(np.abs(pi1024 - uniform1024))

print(f"Stationary distribution:")
print(f"  L1 deviation from uniform: {l1_1024:.6f}")
print(f"  Max deviation from uniform: {max_dev_1024:.6f} ({max_dev_1024/uniform1024*100:.2f}%)")
print()

bset_weight_1024 = sum(pi1024[idx1024[r]] for r in BSet_1024)
print(f"BSet weight: {bset_weight_1024:.6f} (uniform: {len(BSet_1024)/N1024:.6f})")
print()

# Ergodic avg k0
avg_k_1024 = sum(pi1024[i] * v2(r+1) for i, r in enumerate(odd_res_1024))
print(f"Ergodic avg k0 (using v2(r+1)): {avg_k_1024:.6f}")
threshold = 2 * math.log(2) / math.log(1.5)
print(f"Threshold: {threshold:.6f}, Gap: {threshold - avg_k_1024:.6f}")
print()

# Summary
print("=" * 70)
print("EXPANDER CONJECTURE — THREE DATA POINTS")
print("=" * 70)
print()
print(f"{'Modulus':<12} {'States':<8} {'Spectral Gap':<16} {'Max Dev Uniform':<20} {'E[k0] approx'}")
print("-" * 70)
print(f"{'Mod-256':<12} {'128':<8} {'0.938189':<16} {'2.30%':<20} {'1.992'}")
print(f"{'Mod-512':<12} {'256':<8} {'0.920260':<16} {'2.01%':<20} {'1.995'}")
print(f"{'Mod-1024':<12} {N1024:<8} {gap1024:<16.6f} {max_dev_1024/uniform1024*100:<20.2f}% {avg_k_1024:.6f}")
print()

if gap1024 >= 0.90:
    print("EXPANDER CONJECTURE CONFIRMED FOR MOD-1024: gap >= 0.90")
elif gap1024 >= 0.85:
    print("Gap slightly below 0.90 but still large — conjecture may hold with weaker constant.")
else:
    print("Gap below 0.85 — conjecture may be false or threshold is wrong.")

print()
print("TREND ANALYSIS:")
gaps = [0.938189, 0.920260, gap1024]
if len(gaps) >= 3:
    delta_12 = gaps[1] - gaps[0]
    delta_23 = gaps[2] - gaps[1]
    print(f"  Gap decrease mod-256→512: {delta_12:.6f}")
    print(f"  Gap decrease mod-512→1024: {delta_23:.6f}")
    if abs(delta_23) < abs(delta_12):
        print("  Deceleration: gap is stabilizing. Consistent with gap -> constant > 0.")
    else:
        print("  Acceleration: gap decreasing faster. May approach 0 eventually.")
