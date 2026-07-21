"""
139_n20_phantom.py
==================
Find the N=20 phantom cycle discovered in script 138 (Part 6).
Also extend the search more carefully.
"""
import math

def v2(x):
    if x == 0: return 999
    c = 0
    while x % 2 == 0: x //= 2; c += 1
    return c

def macro_step(n):
    K = v2(n+1); m = (n+1) >> K; x = m*(3**K)-1; l0 = v2(x)
    return x >> l0, K, l0

def macro_step_mod(r, N):
    M = 1 << N
    K = v2(r+1); m = (r+1) >> K; x = m*(3**K)-1; l0 = v2(x)
    return (x >> l0) % M, K, l0

def find_cycles_and_phantoms(N):
    M = 1 << N
    odd_residues = list(range(1, M, 2))
    visited = {}
    cycles = []
    for start in odd_residues:
        if start in visited: continue
        path = []; path_set = {}
        n = start
        while n not in visited and n not in path_set:
            path_set[n] = len(path); path.append(n)
            n, _, _ = macro_step_mod(n, N)
        if n in path_set:
            cycle_start = path_set[n]
            cycle = path[cycle_start:]
            cycle_key = tuple(sorted(cycle))
            if cycle_key not in [tuple(sorted(c)) for c in cycles]:
                cycles.append(cycle)
        for node in path:
            visited[node] = True
    phantoms = [c for c in cycles if 1 not in c]
    return phantoms

CRITICAL_RATIO = math.log(4/3) / math.log(2)

print("=" * 70)
print("PART 1: FINDING THE N=20 PHANTOM")
print("=" * 70)
print()
print("Searching for phantom cycles at N=20 (mod 2^20 = 1048576)...")
print("This has 524288 states — using efficient rho algorithm...")

N = 20
M = 1 << N
# Efficient cycle detection without storing all states
# Use Floyd's / functional graph approach

# Actually, we can use a smarter approach: find the longest-orbit element
# and check if it cycles. But with 512K states, let's just do it.
import sys

# Find phantom cycles
phantoms = find_cycles_and_phantoms(N)

print(f"Found {len(phantoms)} phantom cycle(s) at N=20:")
for i, cycle in enumerate(phantoms):
    print(f"\nPhantom {i+1}: length {len(cycle)}")
    print(f"  Elements: {sorted(cycle)[:20]}{'...' if len(cycle)>20 else ''}")

    # Compute K/l0 ratio
    K_vals = []; l0_vals = []
    for r in cycle:
        _, K, l0 = macro_step_mod(r, N)
        K_vals.append(K); l0_vals.append(l0)

    sum_K = sum(K_vals); sum_l0 = sum(l0_vals)
    ratio = sum_l0 / sum_K
    log_gain = sum(k*math.log(3) - (k+l)*math.log(2) for k,l in zip(K_vals,l0_vals))
    net_gain = math.exp(log_gain)

    print(f"  sum(K)={sum_K}, sum(l0)={sum_l0}, ratio={ratio:.5f}")
    print(f"  Critical ratio = {CRITICAL_RATIO:.5f}, imbalance = {ratio-CRITICAL_RATIO:+.5f}")
    print(f"  Net gain per cycle = {net_gain:.4f} ({'expanding' if net_gain > 1 else 'contracting'})")

    # Show dissolution at N=21
    print(f"\n  Dissolution at N+1=21:")
    for r in sorted(cycle)[:10]:
        n_out_exact, K, l0 = macro_step(r)
        n_out_N = n_out_exact % M
        n_out_N1 = n_out_exact % (M << 1)
        in_N = n_out_N in set(cycle)
        in_N1 = n_out_N1 in set(cycle)
        print(f"    f({r}) = {n_out_exact}: mod 2^20 = {n_out_N} (in_cycle={in_N}), "
              f"mod 2^21 = {n_out_N1} (still_in={in_N1})")

    # Smallest element and its real orbit length
    min_r = min(cycle)
    n = min_r; t = 0; max_n = min_r
    while n > 1 and t < 10000:
        n_out, _, _ = macro_step(n)
        n = n_out; t += 1
        if n > max_n: max_n = n
    print(f"\n  Smallest element: {min_r}")
    print(f"  Real orbit length: {t} steps")
    print(f"  Max value in orbit: {max_n} ({max_n.bit_length()} bits)")
    print(f"  Start was: {min_r.bit_length()} bits")

print()
print("=" * 70)
print("PART 2: SYSTEMATIC PHANTOM SURVEY N=11-21")
print("=" * 70)
print()
print("(Checking each N efficiently)")
print()
print(f"{'N':>4} {'#states':>10} {'#phantoms':>12} {'phantom lengths'}")
print("-" * 60)

for N in range(11, 22):
    M = 1 << N
    n_states = M // 2
    if n_states > 1100000:
        print(f"{N:>4} {n_states:>10} {'(too large)':>12}")
        continue

    phantoms = find_cycles_and_phantoms(N)
    lengths = [len(c) for c in phantoms]
    print(f"{N:>4} {n_states:>10} {len(phantoms):>12}  {sorted(lengths) if phantoms else '[]'}")

print()
print("=" * 70)
print("PART 3: K/l0 RATIO ACROSS ALL FOUND PHANTOMS")
print("=" * 70)
print()

# Collect all phantoms across all N
all_phantom_data = {}
# Previously found:
all_phantom_data[7] = [[47, 91, 103, 121]]
all_phantom_data[8] = [[71, 91, 103, 121, 175, 189]]
all_phantom_data[9] = [[91, 95, 103, 167, 175, 253, 283, 319, 399, 445]]
all_phantom_data[10] = [[703, 937]]

print(f"{'N':>4} {'len':>6} {'sum_K':>8} {'sum_l0':>8} {'ratio':>8} {'crit':>8} {'imbal':>8} {'gain':>10}")
print("-" * 70)

for N in sorted(all_phantom_data.keys()):
    for cycle in all_phantom_data[N]:
        K_vals = []; l0_vals = []
        for r in cycle:
            _, K, l0 = macro_step_mod(r, N)
            K_vals.append(K); l0_vals.append(l0)
        sum_K = sum(K_vals); sum_l0 = sum(l0_vals)
        ratio = sum_l0/sum_K
        log_gain = sum(k*math.log(3)-(k+l)*math.log(2) for k,l in zip(K_vals,l0_vals))
        print(f"{N:>4} {len(cycle):>6} {sum_K:>8} {sum_l0:>8} {ratio:>8.5f} "
              f"{CRITICAL_RATIO:>8.5f} {ratio-CRITICAL_RATIO:>+8.5f} {math.exp(log_gain):>10.4f}")

print()
print("Hypothesis: All phantom cycles have gain != 1 (ruling them out as genuine cycles)")
print("If any phantom had gain = 1 exactly, it would BE a genuine Collatz cycle.")
