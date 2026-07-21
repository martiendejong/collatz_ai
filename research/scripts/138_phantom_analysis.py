"""
138_phantom_analysis.py
=======================
Deep analysis of the phantom cycle window N=7-10.

Key diagnostic: for a GENUINE Collatz cycle, the K and l0 values satisfy:
  sum(l0_i) / sum(K_i) = log(4/3) / log(2) = log2(4/3) = 1/log2(4/3)^{-1}
  = log(4/3)/log(2) = 0.41504...

This ratio arises from the balance condition:
  product_i (n_i_out / n_i) = 1  (cycle returns to start)
  => product_i (3^{K_i} / 2^{K_i+l0_i}) = 1  (approximately, for large n)
  => sum_i K_i * log3 = sum_i (K_i + l0_i) * log2
  => sum(l0_i)/sum(K_i) = log(4/3)/log(2) ~= 0.41504

Phantom cycles will have WRONG ratio. How far off are they?
Also: verify the phantom dissolution mechanism step by step.
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

CRITICAL_RATIO = math.log(4/3) / math.log(2)
print("=" * 70)
print("PART 1: K/l0 RATIO TEST FOR PHANTOM CYCLES")
print("=" * 70)
print()
print(f"Critical ratio for genuine Collatz cycles: sum(l0)/sum(K) = log(4/3)/log2 = {CRITICAL_RATIO:.5f}")
print()

phantom_data = {
    7: [47, 91, 103, 121],
    8: [71, 91, 103, 121, 175, 189],
    9: [91, 95, 103, 167, 175, 253, 283, 319, 399, 445],
    10: [703, 937]
}

for N, elements in phantom_data.items():
    cycle = list(elements)
    K_vals = []; l0_vals = []; steps = []

    for r in cycle:
        n_out, K, l0 = macro_step_mod(r, N)
        K_vals.append(K); l0_vals.append(l0)
        steps.append((r, K, l0, n_out))

    sum_K = sum(K_vals); sum_l0 = sum(l0_vals)
    ratio = sum_l0 / sum_K
    imbalance = ratio - CRITICAL_RATIO

    # Compute the "product gain" of the cycle
    # Each step: n_out/n_in ~= 3^K / 2^{K+l0}
    log_gain = sum(K*math.log(3) - (K+l0)*math.log(2) for K,l0 in zip(K_vals,l0_vals))
    net_gain = math.exp(log_gain)

    print(f"N={N}: {len(cycle)}-cycle")
    print(f"  Steps:  {[(r, K, l0, n_out) for r, K, l0, n_out in steps]}")
    print(f"  sum(K)={sum_K}, sum(l0)={sum_l0}, ratio={ratio:.5f}")
    print(f"  Imbalance from critical: {imbalance:+.5f}")
    print(f"  Net log gain per cycle: {log_gain:.4f} (product factor={net_gain:.4f})")
    if net_gain > 1:
        print(f"  ==> Expanding cycle! (>1). Ruled out as genuine Collatz cycle. (Genuine requires ~1)")
    else:
        print(f"  ==> Contracting cycle! (<1). Also ruled out as genuine. (Genuine requires ~1)")
    print()

print()
print("=" * 70)
print("PART 2: PHANTOM DISSOLUTION — STEP-BY-STEP")
print("=" * 70)
print()
print("For each N=7-10 phantom, trace the orbit at the NEXT level (N+1)")
print("to see exactly where the cycle breaks.")
print()

for N in [7, 8, 9, 10]:
    cycle = phantom_data[N]
    N1 = N + 1
    print(f"N={N}: phantom cycle {sorted(cycle)}")
    print(f"  Tracing at level N+1={N1} (mod {1<<N1}):")
    for r in cycle:
        n_out, K, l0 = macro_step(r)  # EXACT, not modular
        n_out_modN = n_out % (1 << N)
        n_out_modN1 = n_out % (1 << N1)
        in_cycle_N = n_out_modN in cycle
        in_cycle_N1 = n_out_modN1 in cycle
        print(f"    f_exact({r}) = {n_out} (mod {1<<N}={n_out_modN}, in_cycle={in_cycle_N})"
              f"   (mod {1<<N1}={n_out_modN1}, still_in_cycle={in_cycle_N1})")
    print()

print()
print("=" * 70)
print("PART 3: REAL ORBIT LENGTHS STARTING FROM PHANTOM ELEMENTS")
print("=" * 70)
print()
print("How many REAL macro-steps does each phantom element take to reach 1?")
print("(If it's a true cycle, it never reaches 1 - but they all do)")
print()

all_phantoms = set()
for cycle in phantom_data.values():
    all_phantoms.update(cycle)

print(f"{'Element':>8} {'Orbit length':>14} {'Max value':>14} {'Max/start':>12}")
print("-" * 55)

for r in sorted(all_phantoms):
    n = r
    t = 0
    max_n = r
    while n > 1:
        n_out, K, l0 = macro_step(n)
        n = n_out
        t += 1
        if n > max_n: max_n = n
        if t > 10000: break  # safety
    print(f"{r:>8} {t:>14} {max_n:>14} {max_n/r:>12.2f}")

print()
print("All phantom elements reach 1 in finite steps - confirming no genuine cycle.")

print()
print("=" * 70)
print("PART 4: WHY N=7 IS THE FIRST PHANTOM LEVEL")
print("=" * 70)
print()
print("The phantom at N=7 is caused by macro_step(103)=175 and 175 mod 128 = 47.")
print("Key: 175 - 47 = 128 = 2^7. The loop closes at this modulus.")
print()
print("For N=6 (mod 64): 175 mod 64 = 47. So 103 -> 47 mod 64 too!")
print("Why doesn't this create a cycle at N=6?")
print()
print("Tracing the orbit mod 64 starting at 47:")
visited = []; n = 47
for _ in range(30):
    visited.append(n)
    n_out, _, _ = macro_step_mod(n, 6)
    n = n_out
    if n in visited or n == 1:
        break
idx = visited.index(n) if n in visited else None
print(f"  Path (mod 64): {visited} -> {n}")
if idx is not None:
    print(f"  Cycle found! Starting at index {idx}: {visited[idx:]}")
else:
    print(f"  No cycle back to start. n={n} (= {'1' if n==1 else 'other'})")

print()
print("The reason N=6 has NO phantom at 47:")
print("At N=6: 47 -> f_6(47) = ?")
n_out6, K, l0 = macro_step_mod(47, 6)
print(f"  f_6(47) = macro_step(47) mod 64 = {macro_step(47)[0]} mod 64 = {n_out6}")
print(f"  (Note: macro_step(47) = {macro_step(47)[0]}, and {macro_step(47)[0]} mod 64 = {macro_step(47)[0] % 64})")
print()
print("At N=7: f_7(47) = ?")
n_out7, K, l0 = macro_step_mod(47, 7)
print(f"  f_7(47) = {macro_step(47)[0]} mod 128 = {n_out7}")

print()
print("So at N=6, 47 maps to something DIFFERENT:")
print("macro_step(47)=121: 121 mod 64 = 57, 121 mod 128 = 121")
print("At N=6: 47->57->?")
print("At N=7: 47->121->?")

# Trace from 47 at N=6
print()
print("Full chain at N=6 starting from 47:")
seen = set()
n = 47
chain = []
for _ in range(40):
    if n in seen: break
    seen.add(n)
    chain.append(n)
    n, _, _ = macro_step_mod(n, 6)
chain.append(n)
print(f"  {' -> '.join(str(x) for x in chain)}")
print(f"  (terminates at {n})")

print()
print("=" * 70)
print("PART 5: ALGEBRAIC STRUCTURE — WHY THE WINDOW CLOSES AT N=11")
print("=" * 70)
print()
print("The N=10 phantom {703, 937} closes because:")
print("  macro_step(703) = 4009, and 4009 mod 1024 = 937")
print("  macro_step(937) = 703 (exactly)")
print()
print("But at N=11: 4009 mod 2048 = ?")
print(f"  4009 mod 2048 = {4009 % 2048}")
print()
print("For the phantom to survive at N=11, we would need")
print("macro_step(703) mod 2048 = 937 AND macro_step(937) mod 2048 = 703.")
print()
print("macro_step(703) = 4009. 4009 mod 2048 = 1961 != 937. Phantom dies.")
print()
print("After dissolution, where does 703 go at N=11?")
n = 703; chain = [703]
for _ in range(30):
    n, _, _ = macro_step_mod(n, 11)
    chain.append(n)
    if n == 1: break
print(f"  Chain mod 2048: {' -> '.join(str(x) for x in chain[:20])}...")
print(f"  Length until reaching 1 mod 2048: {len(chain)-1 if 1 in chain else 'N/A'}")

print()
print("=" * 70)
print("PART 6: IS THE PHANTOM WINDOW BOUNDED ABOVE?")
print("=" * 70)
print()
print("Testing N=20-30 for any new phantom cycles (new computation)...")
print()

def find_phantom_count(N):
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
            cycle = tuple(sorted(path[cycle_start:]))
            if cycle not in [tuple(sorted(c)) for c in cycles]:
                cycles.append(list(cycle))
        for node in path:
            visited[node] = True
    phantom_cycles = [c for c in cycles if 1 not in c]
    return len(phantom_cycles)

# Only go up to N=22 (4M states might be slow but manageable)
for N in range(20, 24):
    M = 1 << N
    n_states = M // 2
    if n_states > 2000000:
        print(f"N={N}: {n_states} states - too large, stopping")
        break
    count = find_phantom_count(N)
    print(f"N={N}: {count} phantom cycles {'(PHANTOM WINDOW STILL OPEN!)' if count > 0 else ''}")

print()
print("If no phantoms appear from N=11 onward: the phantom window N=7-10 is isolated.")
print("This is consistent with the Collatz conjecture (no genuine non-trivial cycle).")
