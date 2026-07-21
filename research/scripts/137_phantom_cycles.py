"""
137_phantom_cycles.py
=====================
The Collatz macro-step functional graph mod 2^N.

Script 136 revealed a shocking pattern:
- For N=4,5,6: only 1 eigenvalue of magnitude 1 (just the true fixed point n=1)
- For N=7,8,9,10: MULTIPLE eigenvalues of magnitude 1 (phantom cycles appear!)
- For N=11,12,13: back to 1 eigenvalue (phantom cycles disappear again)

This means: for some N, the modular map f: (odd mod 2^N) -> (odd mod 2^N)
has PHANTOM CYCLES — cycles that don't correspond to real Collatz cycles.

Key question: What are these phantom cycles?
- Find them explicitly for N where they exist
- Determine if they correspond to hypothetical Collatz cycles that have been
  ruled out by computation (they must, since Collatz verified to ~10^20)
- Understand the algebraic structure creating them

The Collatz conjecture is equivalent to: for ALL N, the ONLY cycle
in the true Collatz map is {1,2,4,8,...}.
But in the MODULAR map mod 2^N, phantom cycles can exist (false positives).
"""
import math

def v2(x):
    if x == 0: return 999
    c = 0
    while x % 2 == 0: x //= 2; c += 1
    return c

def macro_step_exact(n):
    """Exact macro step (not modular)."""
    K = v2(n+1)
    m = (n+1) >> K
    x = m * (3**K) - 1
    l0 = v2(x)
    return x >> l0, K, l0

def macro_step_mod(r, N):
    """Compute macro_step(r) mod 2^N for odd r.
    Only exact if v2(m*3^K - 1) < N - K. Otherwise: approximation."""
    M = 1 << N
    K = v2(r+1)
    m = (r+1) >> K
    x = m * (3**K) - 1
    l0 = v2(x)
    n_out = x >> l0
    return n_out % M

def find_cycles_mod_N(N):
    """Find all cycles in the functional graph of macro_step mod 2^N."""
    M = 1 << N
    odd_residues = list(range(1, M, 2))

    # Follow each starting point until revisiting
    visited = {}  # r -> step_count when first visited
    in_cycle = set()
    cycles = []

    for start in odd_residues:
        if start in visited:
            continue

        # Trace the orbit until we hit a visited node or revisit
        path = []
        path_set = {}
        n = start

        while n not in visited and n not in path_set:
            path_set[n] = len(path)
            path.append(n)
            n = macro_step_mod(n, N)

        if n in path_set:
            # Found a new cycle
            cycle_start_idx = path_set[n]
            cycle = path[cycle_start_idx:]
            cycle_tuple = tuple(sorted(cycle))
            if cycle_tuple not in [tuple(sorted(c)) for c in cycles]:
                cycles.append(cycle)
            for node in cycle:
                in_cycle.add(node)
            # Mark all path nodes as visited
            for node in path:
                visited[node] = True
        elif n in visited:
            # Merged into already-visited component
            for node in path:
                visited[node] = True

    return cycles

print("=" * 70)
print("PART 1: CYCLES IN FUNCTIONAL GRAPH mod 2^N")
print("=" * 70)
print()
print(f"{'N':>4} {'#cycles':>8} {'cycle lengths':>30} {'cycle representatives'}")
print("-" * 80)

all_cycles_by_N = {}

for N in range(3, 20):
    cycles = find_cycles_mod_N(N)
    all_cycles_by_N[N] = cycles
    lengths = sorted([len(c) for c in cycles])
    reps = [min(c) for c in cycles]  # smallest element of each cycle

    cycles_str = str(lengths)
    reps_str = str(sorted(reps))

    if len(cycles) > 1:
        print(f"{N:>4} {len(cycles):>8} {cycles_str:>30}  ** PHANTOM CYCLES **")
    else:
        print(f"{N:>4} {len(cycles):>8} {cycles_str:>30}")

    if len(cycles) > 1:
        print(f"     Cycle reps: {reps_str}")

print()
print("=" * 70)
print("PART 2: EXPLICIT PHANTOM CYCLES")
print("=" * 70)
print()

for N in range(3, 20):
    cycles = all_cycles_by_N[N]
    if len(cycles) > 1:
        print(f"N={N}: {len(cycles)} cycles found")
        for i, cycle in enumerate(cycles):
            is_trivial = 1 in cycle
            print(f"  Cycle {i+1} (length {len(cycle)}): {sorted(cycle)[:10]}... {'<-- trivial (contains 1)' if is_trivial else '<-- PHANTOM'}")
            if not is_trivial:
                # Verify it's actually a cycle
                c = sorted(cycle)
                for r in c:
                    out = macro_step_mod(r, N)
                    print(f"    f({r}) = {out} (in cycle: {out in set(cycle)})")
        print()

print()
print("=" * 70)
print("PART 3: ALGEBRAIC ANALYSIS OF PHANTOM CYCLES")
print("=" * 70)
print()
print("A phantom cycle element r satisfies: f^k(r) = r mod 2^N for some k.")
print("For a 1-cycle (fixed point mod 2^N): f(r) = r mod 2^N")
print()
print("For 2^N-modular fixed points (k=1):")
print("  macro_step(n) ≡ n mod 2^N")
print("  => (m * 3^K - 1) / 2^l0 ≡ n mod 2^N")
print("  => m * 3^K - 1 ≡ n * 2^l0 mod 2^{N+l0}")
print("  => m * 3^K ≡ 1 + n * 2^l0 mod 2^{N+l0}")
print("  where m = (n+1)/2^K, so n = m*2^K - 1")
print("  => m * 3^K ≡ 1 + (m*2^K - 1) * 2^l0 mod 2^{N+l0}")
print("  => m * 3^K ≡ 1 + m * 2^{K+l0} - 2^l0 mod 2^{N+l0}")
print("  => m * (3^K - 2^{K+l0}) ≡ 1 - 2^l0 mod 2^{N+l0}")
print("  => m ≡ (1 - 2^l0) / (3^K - 2^{K+l0}) mod 2^{N+l0} / gcd(3^K - 2^{K+l0}, 2^{N+l0})")
print("  (if 3^K - 2^{K+l0} < 0, there's no positive solution)")
print()

# Find fixed points explicitly
print("Fixed points of macro_step (solutions to f(n) = n):")
print("  K=1, l0=1: m = (2-1)/(3-4) = 1/(-1) -- no positive solution")
print("  K=1, l0=1 (fixing sign): 3^1 - 2^{1+1} = 3 - 4 = -1")
print("  m(3^K - 2^{K+l0}) = 1 - 2^l0")
print("  m(-1) = -1 => m = 1 => n = 2^1 * 1 - 1 = 1 ✓")
print()

# Verify K=1, l0=1 gives n=1
n = 1
K = v2(n+1)  # v2(2) = 1
m = (n+1) >> K  # m = 1
x = m * 3**K - 1  # 3 - 1 = 2
l0 = v2(x)  # v2(2) = 1
n_out = x >> l0  # 2 >> 1 = 1
print(f"  Verify: macro_step(1): K={K}, m={m}, x={x}, l0={l0}, n_out={n_out} ✓")
print()

print("=" * 70)
print("PART 4: TRAJECTORY LENGTHS IN FUNCTIONAL GRAPH")
print("=" * 70)
print()
print("How many steps until each state reaches the cycle?")
print("(This is the 'tail length' in the rho structure)")

def trajectory_to_cycle(r, N, cycle_set, max_steps=100):
    """Number of steps from r to enter the cycle."""
    n = r
    for t in range(max_steps):
        if n in cycle_set:
            return t
        n = macro_step_mod(n, N)
    return max_steps

for N in [6, 8, 10, 12]:
    cycles = all_cycles_by_N[N]
    M = 1 << N
    odd_residues = list(range(1, M, 2))
    cycle_set = set()
    for c in cycles:
        cycle_set.update(c)

    tail_lengths = [trajectory_to_cycle(r, N, cycle_set) for r in odd_residues]
    import statistics
    print(f"N={N}: mean tail = {statistics.mean(tail_lengths):.2f}, max tail = {max(tail_lengths)}, "
          f"fraction in cycle = {len(cycle_set)/len(odd_residues):.4f}")

print()
print("=" * 70)
print("PART 5: 2-ADIC LIFTING OF PHANTOM CYCLES")
print("=" * 70)
print()
print("If r is a phantom cycle element mod 2^N, can it lift to mod 2^{N+1}?")
print("A phantom that survives all liftings would give a true Collatz cycle.")
print("Collatz is verified to n ~ 10^20, so all phantoms must die at some point.")
print()

def find_phantom_elements(N):
    """Find non-trivial cycle elements mod 2^N (elements not in the 1-cycle)."""
    cycles = all_cycles_by_N.get(N, find_cycles_mod_N(N))
    phantom = []
    for c in cycles:
        if 1 not in c:
            phantom.extend(c)
    return set(phantom)

for N in range(3, 18):
    if N not in all_cycles_by_N:
        all_cycles_by_N[N] = find_cycles_mod_N(N)
    phantoms = find_phantom_elements(N)
    if phantoms:
        print(f"N={N}: phantom elements: {sorted(phantoms)[:20]}")

print()
print("=" * 70)
print("PART 6: PATTERN — DO PHANTOM CYCLES FOLLOW A HENSEL-LIKE STRUCTURE?")
print("=" * 70)
print()
print("Hensel's lemma: a root mod p can lift to a root mod p^N.")
print("For Collatz: a phantom cycle mod 2^N might lift to 2^{N+1}, 2^{N+2}, ...")
print("If it lifts to ALL 2^N, it gives a 2-adic cycle (actual Collatz cycle).")
print()

# Check: for each N where phantoms exist, do the phantoms mod 2^N come from
# phantom elements at N-1, N+1?
for N in range(5, 18):
    if N not in all_cycles_by_N:
        all_cycles_by_N[N] = find_cycles_mod_N(N)
    phantoms_N = find_phantom_elements(N)
    if phantoms_N and N+1 <= 17:
        if N+1 not in all_cycles_by_N:
            all_cycles_by_N[N+1] = find_cycles_mod_N(N+1)
        phantoms_N1 = find_phantom_elements(N+1)

        # Check which phantoms mod 2^N lift to phantoms mod 2^{N+1}
        # A phantom r mod 2^N lifts if either r or r+2^N is also a phantom mod 2^{N+1}
        lifts = set()
        for r in phantoms_N:
            if r in phantoms_N1 or (r + (1<<N)) in phantoms_N1:
                lifts.add(r)

        print(f"N={N}: {len(phantoms_N)} phantom elements, {len(lifts)} lift to N={N+1}")
        if phantoms_N:
            lift_fraction = len(lifts) / len(phantoms_N)
            print(f"  Lift fraction: {lift_fraction:.2f}")
