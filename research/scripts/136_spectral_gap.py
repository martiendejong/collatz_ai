"""
136_spectral_gap.py
===================
Spectral gap of the Collatz Markov chain on residues mod 2^N.

The Collatz macro-step induces a Markov chain on odd residues mod 2^N.
From Obs 274 (orbit coupling): 2-adic distance contracts by ~4 per step.
This gives a mixing time T_mix ~ N/4 steps.

Standard spectral theory: T_mix ~ 1/gap where gap = 1 - |lambda_2|.
So: gap ~ 4/N → gap * N → 4 as N → infinity.

Key questions:
1. Does gap ~ C/N (linear decay)?
2. What is the constant C?
3. Is there an exact formula?

Method: For each N, compute the transition matrix P on odd residues mod 2^N,
find all eigenvalues, extract spectral gap = 1 - |lambda_2|.

This is expensive (matrix of size 2^{N-1} x 2^{N-1}), so we go up to N=14.
"""
import numpy as np
from scipy import linalg
import math
from collections import defaultdict

def v2(x):
    if x == 0: return 999
    c = 0
    while x % 2 == 0: x //= 2; c += 1
    return c

def macro_step_mod(n, modulus):
    """Compute macro_step(n) mod modulus for odd n."""
    K = v2(n+1)
    m = (n+1) >> K
    # x = m * 3^K - 1
    x = (m * pow(3, K, modulus * (1 << 100))) - 1  # keep enough precision
    l0 = v2(x)
    n_out = x >> l0
    return n_out % modulus

def build_transition_matrix(N):
    """Build transition matrix for Collatz on odd residues mod 2^N."""
    M = 1 << N  # modulus = 2^N
    # Odd residues mod 2^N: 1, 3, 5, ..., M-1 (total 2^{N-1} = M/2 states)
    odd_residues = list(range(1, M, 2))
    n_states = len(odd_residues)
    state_to_idx = {r: i for i, r in enumerate(odd_residues)}

    # For each state r, the macro_step output n_out depends on n ~ r mod 2^N.
    # The EXACT output depends on n (not just n mod 2^N), but the output mod 2^{N-K}
    # is determined by n mod 2^N (where K = v2(n+1)).
    # However: from n mod 2^N, we know K = v2(n+1 mod 2^N) exactly if n+1 < 2^N.
    # For the modular chain, we need to be careful.

    # Simpler approach: empirical transition probabilities via large orbit sampling.
    # Count transitions: for n ≡ r mod 2^N (uniformly sampled), where does n_out end up?

    # For small N, can enumerate all odd residues mod 2^N directly.
    # macro_step(n) for n ≡ r mod 2^N:
    #   K = v2(n+1). n+1 ≡ r+1 mod 2^N. v2(r+1) = v2(r+1 mod 2^N).
    #   If K < N: K is determined exactly by r. If K >= N: K depends on higher bits.
    # For exact computation: assume K = v2(r+1) if v2(r+1) < N, else K is random.
    # But for the modular chain, we need to marginalize over the unknown higher bits.

    # Method: for each r, enumerate all n in [1, 2^N * S) with n ≡ r mod 2^N,
    # compute macro_step(n) mod 2^N, and estimate transition probabilities.

    # Simplified: use the EXACT modular map when K < N.
    # When K >= N (i.e., v2(r+1) = v2(0+1) = ... hmm, v2(r+1) depends on r specifically)
    # v2(r+1) can be at most N-1 for odd r in [1, 2^N-1] (since r+1 is even, v2(r+1) = v2(r+1))
    # Actually v2(r+1) is well-defined for each r independently.

    # For N up to 14, use direct computation:
    # For each odd r mod 2^N: K = v2(r+1). K is at most N-1 (since r+1 <= 2^N).
    # So K < N always. ✓ (because r <= 2^N-1, r+1 <= 2^N, v2(r+1) <= N but if r+1=2^N then r=2^N-1 which is odd and r+1=2^N so K=N. But r must be < 2^N so r+1 <= 2^N, and r+1=2^N iff r=2^N-1.)
    # For r = 2^N - 1: K = N. This is a special case.

    # Build exact transition matrix (deterministic):
    P = np.zeros((n_states, n_states))

    for i, r in enumerate(odd_residues):
        K = v2(r+1)
        m = (r+1) >> K
        # x = m * 3^K - 1. But we want x mod (2^N * something) for the l0 computation.
        # The issue: l0 = v2(x) and x can depend on higher bits of m.
        # For the MODULAR chain: m = (n+1)/2^K where n ≡ r mod 2^N.
        # m mod 2^{N-K} is determined by r: m ≡ (r+1)/2^K mod 2^{N-K}.
        m_mod = (r+1) >> K  # This is m mod 2^{N-K} (since r+1 < 2^N+1, r+1 >> K < 2^{N-K+1})
        # Actually m = (r+1)/2^K exactly (since v2(r+1) = K and r < 2^N).
        # So m = (r+1) >> K is exact. But x = m * 3^K - 1 is then EXACT (no modular issues).
        x = m * (3**K) - 1
        l0 = v2(x)
        n_out = (x >> l0) % M
        # n_out might be even (if the macro step gave an even output) - shouldn't happen
        if n_out % 2 == 0:
            # This can happen if x >> l0 is even, which would mean we haven't removed all 2s
            # But l0 = v2(x) should give x >> l0 odd. Let's verify.
            pass  # shouldn't happen
        j = state_to_idx.get(n_out, None)
        if j is not None:
            P[i, j] = 1.0
        else:
            # n_out not in our state space (shouldn't happen if mod 2^N is right)
            # Might happen if n_out > 2^N (needs reduction mod 2^N)
            n_out_red = n_out % M
            if n_out_red % 2 == 1:
                j = state_to_idx.get(n_out_red, None)
                if j is not None:
                    P[i, j] = 1.0

    return P, odd_residues

print("=" * 70)
print("PART 1: SPECTRAL GAP AS FUNCTION OF N")
print("=" * 70)
print()
print("For each N, building transition matrix on odd residues mod 2^N...")
print()
print(f"{'N':>5} {'#states':>8} {'|lam_2|':>10} {'gap':>10} {'gap*N':>10} {'T_mix_est':>12}")
print("-" * 60)

gaps = []
N_vals = []

for N in range(4, 16):
    M = 1 << N
    n_states = M // 2
    if n_states > 8000:
        print(f"N={N}: {n_states} states — too large for dense matrix, skipping")
        continue

    P, residues = build_transition_matrix(N)

    # Check if P is stochastic (each row sums to 1)
    row_sums = P.sum(axis=1)
    if not np.allclose(row_sums, 1, atol=0.01):
        print(f"N={N}: WARNING — some rows don't sum to 1: min={row_sums.min():.3f}")

    # Get eigenvalues
    try:
        eigvals = linalg.eigvals(P)
        eigvals_real = np.sort(np.abs(eigvals.real))[::-1]
        lam1 = eigvals_real[0]  # should be ~1
        lam2 = eigvals_real[1]  # second largest
        gap = 1 - lam2
        T_mix = -1/math.log(lam2) if lam2 > 0 and lam2 < 1 else float('inf')
        print(f"{N:>5} {n_states:>8} {lam2:>10.6f} {gap:>10.6f} {gap*N:>10.4f} {T_mix:>12.4f}")
        gaps.append(gap)
        N_vals.append(N)
    except Exception as e:
        print(f"N={N}: Error computing eigenvalues: {e}")

print()
print("=" * 70)
print("PART 2: FITTING gap ~ C/N^alpha")
print("=" * 70)
print()

if len(N_vals) >= 4:
    log_N = np.log(N_vals)
    log_gap = np.log(gaps)
    alpha, log_C = np.polyfit(log_N, log_gap, 1)
    C = math.exp(log_C)
    print(f"Power law fit: gap = {C:.4f} / N^{-alpha:.4f}")
    print(f"(or: gap ~ N^{alpha:.4f})")
    print()
    print("Residuals:")
    for N, g in zip(N_vals, gaps):
        pred = C * N**alpha
        print(f"  N={N:>3}: gap={g:.6f}, pred={pred:.6f}, ratio={g/pred:.4f}")

    # Also fit gap * N^alpha = C
    print(f"\nBest fit alpha = {alpha:.4f}")
    print(f"If alpha = -1 (gap ~ 1/N): gap*N should be constant = {C:.4f}")
    print(f"Observed gap*N values: {[round(g*N,3) for g,N in zip(gaps,N_vals)]}")

print()
print("=" * 70)
print("PART 3: EIGENVALUE SPECTRUM STRUCTURE")
print("=" * 70)
print()
print("Top 10 eigenvalues (by magnitude) for N=10:")

N = 10
P, residues = build_transition_matrix(N)
eigvals = linalg.eigvals(P)
eigvals_sorted = sorted(eigvals, key=lambda x: -abs(x))

print(f"N={N}: {len(residues)} states")
print(f"{'Rank':>6} {'|lam|':>10} {'Re(lam)':>12} {'Im(lam)':>12}")
for i, ev in enumerate(eigvals_sorted[:20]):
    print(f"{i:>6} {abs(ev):>10.6f} {ev.real:>12.6f} {ev.imag:>12.6f}")

print()
print("Gap = 1 - |lam_2|:", 1 - abs(eigvals_sorted[1]))
print()

# Check if eigenvalue 1 has multiplicity 1 (unique stationary distribution)
lam1_count = sum(1 for ev in eigvals_sorted if abs(ev.real - 1) < 0.001 and abs(ev.imag) < 0.001)
print(f"Number of eigenvalues ≈ 1: {lam1_count}")

print()
print("=" * 70)
print("PART 4: STATIONARY DISTRIBUTION vs UNIFORM")
print("=" * 70)
print()
print("Is the stationary distribution of the mod-2^N chain exactly uniform?")

for N in [6, 8, 10]:
    M = 1 << N
    P, residues = build_transition_matrix(N)
    # Find stationary dist (left eigenvector of eigenvalue 1)
    eigvals, leftvecs = linalg.eig(P.T)
    # Find eigenvector for eigenvalue ~1
    idx = np.argmin(abs(eigvals - 1))
    pi = leftvecs[:, idx].real
    pi = pi / pi.sum()
    pi = pi * len(residues)  # normalize to n_states * pi = 1 means uniform

    # Check if uniform
    max_dev = max(abs(pi - 1))
    mean_pi = np.mean(pi)
    std_pi = np.std(pi)
    print(f"N={N}: stationary dist (scaled): mean={mean_pi:.4f}, std={std_pi:.6f}, max_dev={max_dev:.6f}")
    print(f"  {'Uniform' if max_dev < 0.01 else 'NON-UNIFORM'}")

print()
print("=" * 70)
print("PART 5: MIXING TIME PREDICTION vs COUPLING TIME")
print("=" * 70)
print()
print("From Obs 274 (orbit coupling): coupling time T_couple ~ 1.1b")
print("b = N for the mod-2^N chain")
print("Theory: T_mix ~ 1.1*N (from coupling)")
print()
print("From spectral theory: T_mix ~ 1/gap")
print("If gap ~ C/N: T_mix ~ N/C")
print()

if len(N_vals) >= 4:
    # Predicted from coupling: T_mix = 1.1 * N
    # Predicted from spectral: T_mix = 1/gap
    # These should agree if gap ~ 1/(1.1*N)
    C_coupling = 1/1.1  # = 0.909
    print(f"Coupling prediction: gap ~ {C_coupling:.3f}/N")
    print(f"Spectral fit: gap ~ {C:.3f}/N^{-alpha:.3f}")
    print()
    if alpha >= -1.1 and alpha <= -0.9:
        print("Results consistent with gap ~ 1/N (T_mix ~ N).")
    else:
        print(f"Discrepancy: alpha = {alpha:.3f} != -1.0")
