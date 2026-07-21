"""
113_bset_chain_structure.py
============================
DEEP ANALYSIS OF THE EMBEDDED BSet MARKOV CHAIN

1. Compute exact 15x15 transition matrix P_BSet(r->r') empirically.
2. Compute Phi(r) = avg k0 per macro-step in excursion from r.
3. Find dominant transition paths (which elements send to which).
4. Connect coset structure to transition probabilities.
5. Explain why r=103 has highest stationary weight (pi=0.123).
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

BSet = {27,55,63,83,95,103,127,159,169,191,207,223,239,253,255}
BList = sorted(BSet)

# Precompute K and j for each BSet element
def get_K_j(r):
    K = v2(r + 1)
    m_red = (r + 1) >> K
    l0 = v2((3**K) * m_red - 1)
    j = 8 - K - l0
    return K, l0, j

bset_meta = {r: get_K_j(r) for r in BList}

print("=" * 70)
print("PART 1: BSet ELEMENT COSET HIERARCHY")
print("=" * 70)
print()
print(f"{'r':>4}  {'K':>2}  {'l0':>3}  {'j=8-K-l0':>10}  {'P(h=1) exact':>14}")
print("-" * 50)
for r in BList:
    K, l0, j = bset_meta[r]
    if j >= 1:
        n_bset_j = sum(1 for r2 in BSet if v2(r2+1) >= j)
        ph1 = f"{n_bset_j}/{1<<(8-j)}"
    else:
        ph1 = f"~15/128 (j={j})"
    print(f"r={r:3d}  K={K}  l0={l0}  j={j:3d}  P(h=1)={ph1}")

print()

# =====================================================================
# PART 2: EMPIRICAL BSet CHAIN TRANSITION MATRIX
# =====================================================================
print("=" * 70)
print("PART 2: EMPIRICAL BSet TRANSITION MATRIX P_BSet(r->r')")
print("=" * 70)
print()

N_SAMP = 4096  # high accuracy

idx = {r: i for i, r in enumerate(BList)}
P_bset = np.zeros((15, 15))
Phi_r = {}
h_distribution = {}

for i, r_start in enumerate(BList):
    K_start = v2(r_start + 1)
    transitions = np.zeros(15)
    phi_vals = []
    h_vals = []

    for k_iter in range(N_SAMP):
        n = r_start + 256 * k_iter
        if v2(n + 1) != K_start:
            continue
        # Trace until next BSet visit
        n_cur = n
        k0_sum = 0
        h = 0
        while True:
            n_next, k0, l0 = macro_step(n_cur)
            k0_sum += k0
            h += 1
            r_next = n_next % 256
            if r_next in BSet:
                transitions[idx[r_next]] += 1
                phi_vals.append(k0_sum / h)
                h_vals.append(h)
                break
            n_cur = n_next

    total = transitions.sum()
    if total > 0:
        P_bset[i] = transitions / total
    Phi_r[r_start] = sum(phi_vals) / len(phi_vals) if phi_vals else 0
    h_distribution[r_start] = h_vals

print(f"Transition matrix computed ({N_SAMP} samples per BSet element).")
print()

# Print the 15x15 matrix (compact)
print("P_BSet(r_row -> r_col):")
header = "        " + " ".join(f"{r:4d}" for r in BList)
print(header)
print("-" * len(header))
for i, r in enumerate(BList):
    row_str = f"r={r:3d}: " + " ".join(f"{P_bset[i,j]:.3f}" if P_bset[i,j] > 0.01 else "  .  " for j in range(15))
    print(row_str)
print()

# =====================================================================
# PART 3: STATIONARY DISTRIBUTION AND Phi COMPUTATION
# =====================================================================
print("=" * 70)
print("PART 3: STATIONARY DISTRIBUTION AND Phi ERGODIC AVERAGE")
print("=" * 70)
print()

# Power iteration for stationary
pi = np.ones(15) / 15
for _ in range(2000):
    pi_new = pi @ P_bset
    if np.max(np.abs(pi_new - pi)) < 1e-10:
        break
    pi = pi_new

print("Stationary distribution of BSet embedded chain:")
print(f"  {'r':>4}  {'K':>2}  {'j':>3}  {'pi':>8}  {'pi*15':>6}  {'Phi(r)':>8}  {'E[h]':>6}")
print("-" * 55)
for i, r in enumerate(BList):
    K, l0, j = bset_meta[r]
    h_vals = h_distribution[r]
    e_h = sum(h_vals)/len(h_vals) if h_vals else 0
    print(f"  r={r:3d}  K={K}  j={j:>3}  pi={pi[i]:.5f}  {pi[i]*15:5.3f}x  Phi={Phi_r[r]:.4f}  E[h]={e_h:.3f}")

print()
print(f"Ergodic avg Phi = {sum(pi[i]*Phi_r[r] for i,r in enumerate(BList)):.6f}")
print(f"Threshold       = {2*math.log(2)/math.log(1.5):.6f}")
print(f"Gap             = {2*math.log(2)/math.log(1.5) - sum(pi[i]*Phi_r[r] for i,r in enumerate(BList)):.6f}")
print()

# =====================================================================
# PART 4: WHY r=103 IS MOST VISITED — FLOW ANALYSIS
# =====================================================================
print("=" * 70)
print("PART 4: FLOW ANALYSIS — WHY r=103 HAS HIGHEST pi")
print("=" * 70)
print()

# For each element r', compute total flow arriving at r' = sum_r pi(r)*P(r->r')
for i_target, r_target in enumerate(BList):
    flow = sum(pi[j] * P_bset[j, i_target] for j in range(15))
    K_t, l0_t, j_t = bset_meta[r_target]
    print(f"  r={r_target:3d} (K={K_t}, j={j_t}): flow_in = {flow:.5f}, pi = {pi[i_target]:.5f}")

print()
print("Top senders to r=103:")
i_103 = idx[103]
contributions = [(P_bset[j, i_103], BList[j]) for j in range(15)]
contributions.sort(reverse=True)
for prob, r_sender in contributions[:8]:
    K_s, l0_s, j_s = bset_meta[r_sender]
    print(f"  r={r_sender:3d} (K={K_s}, j={j_s}): P(send to 103) = {prob:.5f}, contrib = {pi[idx[r_sender]] * prob:.5f}")

print()

# =====================================================================
# PART 5: COSET PREDICTION VS EMPIRICAL
# =====================================================================
print("=" * 70)
print("PART 5: COSET STRUCTURE PREDICTION VS EMPIRICAL TRANSITIONS")
print("=" * 70)
print()
print("For j>=1 elements: predicted P_BSet(r->r') based on coset formula.")
print("Direct transitions (h=1): uniform over K'>=j elements in BSet.")
print("Indirect transitions (h>1): approximated as uniform over all 15.")
print()

for i, r in enumerate(BList):
    K, l0, j = bset_meta[r]
    if j < 1:
        continue
    n_direct = sum(1 for r2 in BSet if v2(r2+1) >= j)
    ph1 = n_direct / (1 << (8-j))
    ph_gt1 = 1 - ph1

    print(f"r={r:3d} (K={K}, j={j}): P(h=1)={ph1:.4f}, P(h>1)={ph_gt1:.4f}")
    print(f"  Predicted: {1/n_direct:.4f} to each of {n_direct} high-K' elements + {ph_gt1/15:.4f} indirect to all")
    print(f"  Empirical: ", end="")
    top_transitions = sorted([(P_bset[i,jj], BList[jj]) for jj in range(15)], reverse=True)[:5]
    for prob, r2 in top_transitions:
        if prob > 0.01:
            K2 = v2(r2+1)
            print(f"r'={r2}(K={K2}):{prob:.3f} ", end="")
    print()

print()

# =====================================================================
# PART 6: SPECTRAL ANALYSIS OF BSet CHAIN
# =====================================================================
print("=" * 70)
print("PART 6: SPECTRAL ANALYSIS OF THE 15x15 BSet CHAIN")
print("=" * 70)
print()

vals_bset = np.linalg.eigvals(P_bset)
vals_bset_sorted = sorted(vals_bset, key=lambda x: -x.real)

print("Eigenvalues of 15x15 BSet embedded chain:")
for i, v in enumerate(vals_bset_sorted[:10]):
    print(f"  lambda_{i+1} = {v.real:.6f} + {v.imag:.4f}i")
print()
gap_bset = 1 - vals_bset_sorted[1].real
print(f"Spectral gap of BSet chain = {gap_bset:.6f}")
print(f"Second eigenvalue = {vals_bset_sorted[1].real:.6f}")
print(f"Mixing time (BSet visits) ~ {1/gap_bset:.2f}")
print()

# The BSet chain spectral gap tells us how fast BSet visits equilibrate
print("INTERPRETATION:")
print(f"  The BSet embedded chain mixes in ~{1/gap_bset:.1f} BSet visits.")
print(f"  With avg E[h] ~ {sum(sum(h_distribution[r]) for r in BList)/sum(len(h_distribution[r]) for r in BList):.2f} macro-steps per BSet visit,")
mean_h = sum(pi[i]*sum(h_distribution[r])/len(h_distribution[r]) for i,r in enumerate(BList) if h_distribution[r])
print(f"  and stationary E[h] = {mean_h:.2f},")
print(f"  full chain mixes in ~{mean_h/gap_bset:.2f} macro-steps.")
