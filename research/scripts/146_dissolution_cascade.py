"""
146_dissolution_cascade.py
===========================
The phantom staircase is a DISSOLUTION CASCADE.

Each phantom cycle at level N has exactly ONE "dissolution point" -- the
element where macro_step(n) > 2^N (the real output exceeds the modulus).

The dissolution cascade:
  N=7 cycle: dissolves at n=103 -> real output 175 = (47 + 2^7) -> ENTERS N=8/9 cycle
  N=8 cycle: dissolves at n=175 -> real output 445 = (189 + 2^8) -> ENTERS N=9 cycle
  N=9 cycle: dissolves at n=319 -> real output 911 = (399 + 2^9) -> EXIT RAMP -> 23-channel
  N=10 cycle: dissolves at n=703 -> real output 4009 = (937 + 3*2^10) -> long path -> 13-channel

Each level's dissolution leads into the next level's cycle (or to a terminal channel),
creating the phantom STAIRCASE observed at T-8 to T-16.

Key claim: the channel assignment (23 vs 13) is determined by which phantom cycle
dissolution the orbit passes through. This explains the complete disjointness of
the two channels' phantom content.
"""
import math

def v2(x):
    c = 0
    while x % 2 == 0: x //= 2; c += 1
    return c

def macro_step(n):
    K = v2(n+1); m = (n+1) >> K; x = m*(3**K)-1; l0 = v2(x)
    return x >> l0, K, l0

def macro_step_mod(r, N):
    M = 1 << N
    K = v2(r+1); m = (r+1) >> K; x = m*(3**K)-1; l0 = v2(x)
    return (x >> l0) % M

def collatz_to_end(n):
    path = [n]
    while n > 1:
        n, _, _ = macro_step(n)
        path.append(n)
    return path

# Known phantom cycles
PHANTOM_CYCLES = {
    7:  [47, 91, 103, 121],
    8:  [71, 91, 103, 121, 175, 189],
    9:  [91, 95, 103, 167, 175, 253, 283, 319, 399, 445],
    10: [703, 937],
}

print("=" * 70)
print("PART 1: DISSOLUTION POINTS OF EACH PHANTOM CYCLE")
print("=" * 70)
print()
print("For each element p in the phantom cycle at level N,")
print("compute macro_step(p) and check if the result exceeds 2^N.")
print("The dissolution point is the unique element where the result > 2^N.")
print()

for N, cycle in sorted(PHANTOM_CYCLES.items()):
    M = 1 << N
    print(f"N={N} phantom cycle: {cycle}")
    print(f"  Modulus 2^{N} = {M}")
    dissolution_found = []
    for p in cycle:
        n_out, K, l0 = macro_step(p)
        n_out_mod = n_out % M
        exceeds = n_out >= M
        phantom_next = n_out_mod  # what the phantom cycle "expects"
        in_cycle = phantom_next in set(cycle)
        print(f"  p={p:>5}: macro_step -> {n_out:>6} (mod {M}: {n_out_mod:>5}), "
              f"exceeds={str(exceeds):>5}, phantom_next={phantom_next} in_cycle={in_cycle}")
        if exceeds:
            dissolution_found.append((p, n_out, phantom_next))
    print(f"  DISSOLUTION POINTS: {dissolution_found}")
    print()

print()
print("=" * 70)
print("PART 2: DISSOLUTION CASCADE -- WHERE DOES EACH DISSOLUTION LEAD?")
print("=" * 70)
print()
print("Following the real orbit from each dissolution point to see which")
print("phantom cycle it enters next (or which terminal channel it reaches).")
print()

ALL_PHANTOM_FLAT = {}
for N, cycle in PHANTOM_CYCLES.items():
    for v in cycle:
        if v not in ALL_PHANTOM_FLAT:
            ALL_PHANTOM_FLAT[v] = []
        ALL_PHANTOM_FLAT[v].append(N)

# For N=7: dissolution at 103 -> 175
# For N=8: dissolution at 175 -> 445
# For N=9: dissolution at 319 -> 911
# For N=10: dissolution at 703 -> 4009

for N, cycle in sorted(PHANTOM_CYCLES.items()):
    M = 1 << N
    for p in cycle:
        n_out, K, l0 = macro_step(p)
        if n_out < M:
            continue  # not the dissolution point
        # Found dissolution at p -> n_out
        n_out_mod = n_out % M  # what phantom cycle shows
        print(f"N={N}: dissolution at p={p}")
        print(f"  Real:    macro_step({p}) = {n_out}")
        print(f"  Phantom: {p} -> {n_out_mod} (mod 2^{N})")
        print(f"  Difference: {n_out} - {n_out_mod} = {n_out - n_out_mod} = {n_out - n_out_mod // M}*2^{N}+{n_out_mod}")
        print()

        # Follow real orbit from n_out to see where it goes
        path = collatz_to_end(n_out)
        T = len(path) - 1
        print(f"  Real orbit from {n_out}: {T} steps to 1")
        print(f"  T-2 value (2 steps before 1): {path[-3] if len(path) >= 3 else 'N/A'}")
        channel = "23" if (len(path) >= 3 and path[-3] == 23) else \
                  ("13" if (len(path) >= 3 and path[-3] == 13) else "other")
        print(f"  Channel: {channel}-channel")
        print()

        # Show where orbit enters a phantom cycle
        print(f"  First phantom elements entered:")
        for i, v in enumerate(path):
            if v in ALL_PHANTOM_FLAT:
                phantom_Ns = ALL_PHANTOM_FLAT[v]
                print(f"    step {i}: n={v} (phantom N={','.join(str(n) for n in phantom_Ns)}) -> "
                      f"T-{T-i} steps from 1")
                if len([x for x in path[i:] if x in ALL_PHANTOM_FLAT]) < 2:
                    break
        print()

print()
print("=" * 70)
print("PART 3: ALGEBRAIC STRUCTURE OF THE CASCADE")
print("=" * 70)
print()
print("Showing that each dissolution is: real_output = phantom_next + c * 2^N")
print("where c is the 'carry' of the modular overflow.")
print()

for N, cycle in sorted(PHANTOM_CYCLES.items()):
    M = 1 << N
    for p in cycle:
        n_out, K, l0 = macro_step(p)
        if n_out < M:
            continue
        n_out_mod = n_out % M
        c = (n_out - n_out_mod) // M
        print(f"N={N}: p={p}")
        print(f"  macro_step({p}) = {n_out}")
        print(f"  = {n_out_mod} + {c} x 2^{N}")
        print(f"  Phantom next: {n_out_mod}")
        print(f"  Carry c = {c}")
        print()

print()
print("=" * 70)
print("PART 4: CHAIN OF DISSOLUTIONS (THE STAIRCASE MECHANISM)")
print("=" * 70)
print()
print("The staircase is formed by the chain of dissolutions:")
print("  N=7: dissolves at 103 -> 175 [enters N=8/9 cycle, continuing the chain]")
print("  N=8: dissolves at 175 -> 445 [enters N=9 cycle]")
print("  N=9: dissolves at 319 -> 911 [exits phantom zone, begins 23-channel exit ramp]")
print("  N=10: dissolves at 703 -> 4009 [exits phantom zone, begins 13-channel long path]")
print()

# Verify the chain: dissolution output enters next cycle
print("Verification:")
print()
print("N=7: 103 -> 175. Is 175 in N=8 or N=9 cycle? ", 175 in PHANTOM_CYCLES[8], "(N=8) /", 175 in PHANTOM_CYCLES[9], "(N=9)")
print("N=8: 175 -> 445. Is 445 in N=9 cycle? ", 445 in PHANTOM_CYCLES[9])
print("N=9: 319 -> 911. Is 911 in any phantom cycle? ", 911 in set(v for c in PHANTOM_CYCLES.values() for v in c))
print()
print("N=9 exit: 911 starts the 23-channel exit ramp.")
path_911 = collatz_to_end(911)
print(f"911 orbit ({len(path_911)-1} steps): {path_911}")
print(f"T-2 value: {path_911[-3]} -> channel: {'23' if path_911[-3]==23 else '13'}")
print()
print("N=10 exit: 4009 starts the 13-channel long path.")
path_703 = collatz_to_end(703)
print(f"703 orbit ({len(path_703)-1} steps): {path_703}")
print(f"T-2 value: {path_703[-3]} -> channel: {'23' if path_703[-3]==23 else '13'}")

print()
print("=" * 70)
print("PART 5: WHICH PHANTOM ELEMENTS ARE 'TRANSIT' vs 'DISSOLUTION' NODES?")
print("=" * 70)
print()
print("For each phantom element p, classify as:")
print("  TRANSIT: macro_step(p) is also in a phantom cycle (orbit continues in phantom zone)")
print("  DISSOLUTION: macro_step(p) exits the phantom zone (real output > 2^N_max for all N)")
print()

ALL_PHANTOM = set(v for c in PHANTOM_CYCLES.values() for v in c)

for N, cycle in sorted(PHANTOM_CYCLES.items()):
    M = 1 << N
    print(f"N={N} phantom cycle:")
    for p in sorted(cycle):
        n_out, K, l0 = macro_step(p)
        n_out_in_phantom = n_out in ALL_PHANTOM
        Ns_in = [N2 for N2, c2 in PHANTOM_CYCLES.items() if n_out in c2]
        modN = p % 3
        if n_out_in_phantom:
            status = f"TRANSIT -> {n_out} [phantom N={','.join(str(n) for n in Ns_in)}]"
        elif n_out < M:
            status = f"INTERNAL -> {n_out} [real stays in range but not in phantom set]"
        else:
            status = f"DISSOLUTION -> {n_out} [exits phantom zone]"
        print(f"  p={p:>5} (mod3={modN}): {status}")
    print()

print()
print("=" * 70)
print("PART 6: PROOF THAT THE STAIRCASE IS EXACT")
print("=" * 70)
print()
print("The canonical terminal path:")
print("47 -> 121 -> 91 -> 103 -> 175 -> 445 -> 167 -> 283 -> 319 -> 911 -> 577 -> 433 -> 325 -> 61 -> 23 -> 5 -> 1")
print()
print("Classifying each step:")
canonical = [47, 121, 91, 103, 175, 445, 167, 283, 319, 911, 577, 433, 325, 61, 23, 5, 1]
for i in range(len(canonical)-1):
    p = canonical[i]
    q = canonical[i+1]
    in_N7 = p in PHANTOM_CYCLES[7]
    in_N8 = p in PHANTOM_CYCLES[8]
    in_N9 = p in PHANTOM_CYCLES[9]
    in_N10 = p in PHANTOM_CYCLES[10]
    phantom_levels = [N for N, c in PHANTOM_CYCLES.items() if p in c]
    is_dissolution = [N for N in phantom_levels if macro_step(p)[0] >= (1 << N)]
    q_in_phantom = q in ALL_PHANTOM
    q_Ns = [N for N, c in PHANTOM_CYCLES.items() if q in c]

    if phantom_levels:
        lvl_str = f"phantom N={','.join(str(n) for n in phantom_levels)}"
        if is_dissolution:
            status = f"DISSOLVES at N={','.join(str(n) for n in is_dissolution)}, exits to {q}"
        else:
            status = f"TRANSIT to {q} [phantom N={','.join(str(n) for n in q_Ns)}]"
    else:
        status = f"EXIT RAMP -> {q}"

    print(f"  {p:>5} [{lvl_str if phantom_levels else 'exit ramp':>20}]: {status}")
