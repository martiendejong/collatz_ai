"""
97_r255_selfloop_exact.py
==========================
Deep analysis of the r=255 self-loop — the BEST CYCLE in the BSet graph.

From script 96: T(255,255)=3.7%, E[h]=8.49, k/step=2.5287 → MCM of BSet graph.

These are Collatz trajectories starting at n≡255 mod 256 (k0=8) that
return to r=255 BEFORE hitting any other BSet element.

KEY QUESTIONS:
1. What is the exact mod-256 arithmetic characterizing these paths?
2. What is the h distribution for r=255→r=255 direct returns?
3. How does k/step break down over h?
4. Is there an algebraic reason why r=255 self-loop beats r=255↔r=127 2-cycle?
5. Can we compute an EXACT upper bound on the r=255 self-loop cycle mean?

EXACT COMPUTATION:
For n≡255 mod 256, the m-class has m odd, n=(256m-1). All 256 odd m in [1,511]
give 256 starting points. For each, trace until BSet is hit, check if it's r=255.
"""
import sys, time
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

# =====================================================================
# PART 1: EXACT SMALL-n ANALYSIS OF r=255 SELF-LOOP
# =====================================================================
print("=" * 70)
print("PART 1: EXACT r=255→r=255 SELF-LOOP (small n, exact arithmetic)")
print("=" * 70)
print()
print("For n=255×m×256-1+256k (n≡255 mod 256), k0=v2(n+1)...")
print("Tracing until BSet hit, checking for r=255 return.")
print()

# Use small n for exact computation (n ≡ 255 mod 256)
# Small n: 255, 511, 767, ..., starting from 255+256*k for k=0,1,...
# But we want k0=8 exactly → n+1 = 256m, m odd.

M_BASE_SMALL = 0  # start from very small n
n_values_255 = [256 * m - 1 for m in range(1, 512, 2)]  # n≡255 mod 256 with k0=8 (approximately)
# Actually: n ≡ 255 mod 256 with k0=v2(n+1)=? For n=256m-1: n+1=256m. v2(256m)=8+v2(m).
# For k0=8 exactly: m must be odd. So n=256×(2j+1)-1 for j=0,1,...
# n=255, 767, 1279, ... (step 512)

n_k8 = [256*(2*j+1)-1 for j in range(256)]  # 256 starting points with k0=8

print(f"Using {len(n_k8)} starting points with exact k0=8 (n=255,767,1279,...)")
print(f"Range: n={n_k8[0]} to n={n_k8[-1]}")
print()

MAX_H = 200
h_dist_self = Counter()        # h distribution for r=255 self-loop
h_dist_other = Counter()       # h distribution for r=255→other BSet
dest_dist = Counter()          # destination distribution
k_by_h = defaultdict(list)    # k/step by h for self-loops

self_loop_data = []  # (h, k_sum, k_seq) for each self-loop trajectory
other_data = []      # (h, k_sum, dest) for non-self-loop trajectories

for n in n_k8:
    n_cur = n
    h = 0
    k_sum = 0
    k_seq = []
    dest = None

    while h < MAX_H:
        n_out, k, l = macro_step(n_cur)
        h += 1
        k_sum += k
        k_seq.append(k)
        n_cur = n_out
        if n_cur <= 1:
            dest = 'converged'
            break
        r_out = n_cur % 256
        if r_out in BSet:
            dest = r_out
            break

    if dest == 255:
        h_dist_self[h] += 1
        self_loop_data.append((h, k_sum, k_seq))
    elif dest is not None and dest != 'converged':
        h_dist_other[h] += 1
        dest_dist[dest] += 1
        other_data.append((h, k_sum, dest))

n_self = len(self_loop_data)
n_total = len(n_k8)
print(f"Self-loop (r=255→r=255) count: {n_self}/{n_total} = {100*n_self/n_total:.2f}%")
print(f"  (T(255,255) ≈ {n_self/n_total:.4f})")
print()

if self_loop_data:
    h_vals = [d[0] for d in self_loop_data]
    k_sums = [d[1] for d in self_loop_data]
    avg_h = sum(h_vals) / len(h_vals)
    total_k = sum(k_sums)
    total_h = sum(h_vals)
    avg_kstep = total_k / total_h
    print(f"Self-loop statistics:")
    print(f"  avg h = {avg_h:.3f}")
    print(f"  avg k/step = {avg_kstep:.6f}  (= {total_k}/{total_h})")
    print(f"  max k/step = {max(k/h for h,k,_ in self_loop_data):.4f}")
    print()

    print("h distribution for r=255 self-loop:")
    for h in sorted(h_dist_self):
        cnt = h_dist_self[h]
        paths = [(k, ks) for hs, ks, k in self_loop_data if hs == h]
        avg_k_for_h = sum(ks for _, ks in paths) / sum(hs for hs, _, _ in [(h, ks, k) for hs, ks, k in self_loop_data if hs == h])
        print(f"  h={h:3d}: {cnt:4d} paths  avg k/step={sum(ks for _,ks in paths)/sum([h]*len(paths)):.4f}")

# =====================================================================
# PART 2: k-SEQUENCE PATTERNS IN SELF-LOOP TRAJECTORIES
# =====================================================================
print()
print("=== PART 2: k-SEQUENCE PATTERNS IN r=255 SELF-LOOP ===\n")
print("First step is always k=8 (since k0=8). What follows?")
print()

# Show the h=1 self-loop (if any): n=255, next step lands at r=255
h1_self = [(h, k, ks) for h, ks, k in self_loop_data if h == 1]
if h1_self:
    print(f"h=1 self-loops: {len(h1_self)}")
    for h, k_seq, ks in h1_self[:3]:
        print(f"  k_seq={k_seq}  k_sum={ks}")
else:
    print("h=1 self-loops: 0 (no direct h=1 return to r=255)")

print()

# Show k-sequences for best self-loops (highest k/step)
best_self = sorted(self_loop_data, key=lambda x: -x[1]/x[0])
print("Top 10 self-loop paths by k/step:")
for h, ks, k_seq in best_self[:10]:
    kstep = ks/h
    print(f"  h={h:3d}: k_sum={ks:5d}  k/step={kstep:.4f}  seq={k_seq}")

# =====================================================================
# PART 3: COMPARE r=255 SELF-LOOP vs r=255→r=127 EDGE
# =====================================================================
print()
print("=== PART 3: r=255 SELF-LOOP vs r=255→r=127 EDGE ===\n")
to_127 = [(h, ks, d) for h, ks, d in other_data if d == 127]
if to_127:
    h_127 = [d[0] for d in to_127]
    k_127 = [d[1] for d in to_127]
    avg_kstep_127 = sum(k_127) / sum(h_127)
    print(f"r=255→r=127 edge: {len(to_127)}/{n_total} paths = {100*len(to_127)/n_total:.2f}%")
    print(f"  avg k/step = {avg_kstep_127:.6f}")
    print()
    print(f"r=255→r=255 self: {n_self}/{n_total} paths = {100*n_self/n_total:.2f}%")
    if self_loop_data:
        print(f"  avg k/step = {avg_kstep:.6f}")
    print()

    # For 2-cycle (255↔127) to beat self-loop:
    # need (k255_127 × h255_127 + k127_255 × h127_255) / (h255_127 + h127_255) > k_selfloop
    # Use small-n exact values for r=127→r=255 from the same computation

    # Compute r=127→r=255 transition
    n_k7 = [128*(2*j+1)-1 for j in range(256)]  # n≡127 mod 256 with k0=7
    to_255_from_127 = []
    for n in n_k7:
        n_cur = n
        h = 0
        k_sum = 0
        dest = None
        while h < MAX_H:
            n_out, k, l = macro_step(n_cur)
            h += 1
            k_sum += k
            n_cur = n_out
            if n_cur <= 1:
                dest = 'converged'
                break
            r_out = n_cur % 256
            if r_out in BSet:
                dest = r_out
                break
        if dest == 255:
            to_255_from_127.append((h, k_sum))

    if to_255_from_127:
        h_127_255 = [d[0] for d in to_255_from_127]
        k_127_255 = [d[1] for d in to_255_from_127]
        avg_kstep_127_255 = sum(k_127_255) / sum(h_127_255)
        avg_h_127_255 = sum(h_127_255) / len(h_127_255)
        print(f"r=127→r=255 edge (exact): {len(to_255_from_127)}/256 paths")
        print(f"  avg k/step = {avg_kstep_127_255:.6f}")
        print(f"  avg h = {avg_h_127_255:.2f}")
        print()

        # 2-cycle mean
        k255_127 = avg_kstep_127  # avg k/step for r=255→r=127
        h255_127 = sum(h_127) / len(h_127)
        k127_255 = avg_kstep_127_255
        h127_255 = avg_h_127_255
        two_cycle_mean = (k255_127 * h255_127 + k127_255 * h127_255) / (h255_127 + h127_255)
        print(f"2-cycle mean (r=255↔r=127):")
        print(f"  (k255_127×h255_127 + k127_255×h127_255) / (h255_127 + h127_255)")
        print(f"  = ({k255_127:.4f}×{h255_127:.2f} + {k127_255:.4f}×{h127_255:.2f}) / ({h255_127:.2f} + {h127_255:.2f})")
        print(f"  = {two_cycle_mean:.6f}")
        print()
        if self_loop_data:
            print(f"r=255 self-loop mean: {avg_kstep:.6f}")
            print(f"r=255↔r=127 2-cycle:  {two_cycle_mean:.6f}")
            print(f"Self-loop wins by:    {avg_kstep - two_cycle_mean:.6f}")
    else:
        print("r=127→r=255: 0 direct transitions found in small-n sample")

# =====================================================================
# PART 4: EXACT MCM BOUND FROM SMALL-n COMPUTATION
# =====================================================================
print()
print("=== PART 4: EXACT MCM FROM SMALL-n COMPLETE ENUMERATION ===\n")
print("Using all n≡255 mod 256 with k0=8 in [255, 255+512×256-1]")
print()

if self_loop_data:
    exact_cycle_mean = total_k / total_h
    print(f"r=255 self-loop (small-n exact):")
    print(f"  n_paths = {n_self}")
    print(f"  total k_sum = {total_k}")
    print(f"  total h = {total_h}")
    print(f"  cycle mean = {total_k}/{total_h} = {exact_cycle_mean:.8f}")
    print()
    print(f"  D_hard_kern threshold: 3.419")
    print(f"  Gap: {3.419 - exact_cycle_mean:.6f}")

# =====================================================================
# PART 5: MOD-256×MOD-256 STRUCTURE OF SELF-LOOPS
# =====================================================================
print()
print("=== PART 5: MOD-256 STRUCTURE OF SELF-LOOP TRAJECTORIES ===\n")
print("For each self-loop, show the sequence of mod-256 residues visited")
print("(the 'fingerprint' of each path)")
print()

if self_loop_data:
    # Re-trace self-loops to get mod-256 sequence
    print("Top 5 self-loops by k/step (mod-256 residue sequences):")
    for h, ks, k_seq in best_self[:5]:
        n = n_k8[0]  # need to find which n gives this path
        # Find the actual starting n for this (h, ks) pair
        for n_start in n_k8:
            n_cur = n_start
            h_trace = 0
            k_sum_trace = 0
            mod_seq = [255]  # starting residue
            while h_trace < h:
                n_out, k, l = macro_step(n_cur)
                h_trace += 1
                k_sum_trace += k
                n_cur = n_out
                r_out = n_cur % 256
                mod_seq.append(r_out)
                if r_out in BSet:
                    break
            if h_trace == h and k_sum_trace == ks:
                kstep = ks/h
                print(f"  h={h}, k/step={kstep:.4f}: mod-256 seq = {mod_seq}")
                break

# =====================================================================
# PART 6: SUMMARY
# =====================================================================
print()
print("=" * 70)
print("SUMMARY: r=255 SELF-LOOP IS THE BEST CYCLE")
print("=" * 70)
print()
if self_loop_data:
    print(f"r=255 self-loop cycle mean (small-n exact): {exact_cycle_mean:.6f}")
print(f"Best single edge (r=255→r=127):             2.8013 (large-n, script 96)")
two_cycle_str = f"{two_cycle_mean:.6f}" if to_127 and to_255_from_127 else 'N/A'
print(f"r=255↔r=127 2-cycle mean:                   {two_cycle_str}")
print()
print(f"D_hard_kern threshold: 3.419")
if self_loop_data:
    print(f"Gap (self-loop):       {3.419 - exact_cycle_mean:.6f}")
print()
print("ALGEBRAIC STRUCTURE: The r=255 self-loop trajectories satisfy:")
print("  n ≡ 255 mod 256, k0=8 for first step")
print("  All subsequent steps avoid BSet mod-256 residues")
print("  Final step lands at n' ≡ 255 mod 256 (k0=8)")
print("  The intermediate path through non-BSet territory has high k/step")
print("  because r=255 (k0=8) naturally produces large k values in its excursion")
