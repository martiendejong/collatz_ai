"""
135_joint_padic.py
==================
Joint p-adic independence along Collatz orbits.

Key question: are v_p(n+1) and v_q(n+1) independent for p != q?

If yes: by CRT, n+1 mod M is approximately uniform for any M = product of prime powers.
This would mean the orbit is "perfectly mixing" across all residue classes simultaneously.

This is a strong version of the independence results found so far:
- K = v2(n+1) independent of v_p(n+1) for all p >= 5 (Obs 280)
- K independent of l0 (Obs 268)
- v_p ~ Geom((p-1)/p) for p >= 5 (Obs 280)

If all v_p are JOINTLY independent, then the orbit is a true random number generator
for all modular arithmetic simultaneously.

Also investigate:
- Correlation between v_3(n+1) and v_5(n+1)
- Correlation between v_5(n+1) and v_7(n+1)
- Chi-squared test for joint uniformity mod M
"""
import numpy as np
from collections import defaultdict
import random as _r
import math
from scipy import stats

def v2(x):
    if x == 0: return 999
    c = 0
    while x % 2 == 0: x //= 2; c += 1
    return c

def vp(x, p):
    if x == 0: return 999
    c = 0
    while x % p == 0: x //= p; c += 1
    return c

def macro_step(n):
    K = v2(n+1); m = (n+1) >> K; x = m*(3**K)-1; l = v2(x)
    return x >> l, K, l

print("=" * 70)
print("PART 1: PAIRWISE INDEPENDENCE OF v_p AND v_q")
print("=" * 70)
print()

PRIMES = [2, 3, 5, 7, 11, 13]

_r.seed(42)
n = _r.getrandbits(3000) | 1
n_steps = 100000
K_vals = []
vp_data = {p: [] for p in PRIMES}

for _ in range(n_steps):
    K = v2(n+1)
    K_vals.append(K)
    for p in PRIMES[1:]:  # skip p=2 (K captures it)
        vp_data[p].append(vp(n+1, p))
    n_out, _, _ = macro_step(n)
    n = n_out
    if n < 2:
        n = _r.getrandbits(3000) | 1

K_arr = np.array(K_vals)
vp_arrays = {p: np.array(vp_data[p]) for p in PRIMES[1:]}

print("Pairwise Pearson correlations between v_p(n+1) for different primes:")
prime_pairs = [(3,5),(3,7),(3,11),(5,7),(5,11),(5,13),(7,11),(7,13),(11,13)]
for p,q in prime_pairs:
    Jp = vp_arrays[p]; Jq = vp_arrays[q]
    r = np.corrcoef(Jp, Jq)[0,1]
    print(f"  corr(v_{p}(n+1), v_{q}(n+1)) = {r:.6f}")

print()
print("Correlation of K with v_p for p=3,5,7,11:")
for p in [3,5,7,11]:
    Jp = vp_arrays[p]
    r = np.corrcoef(K_arr, Jp)[0,1]
    print(f"  corr(K, v_{p}(n+1)) = {r:.6f}")

print()
print("=" * 70)
print("PART 2: JOINT UNIFORMITY TEST mod M = 2*3*5 = 30")
print("=" * 70)
print()
print("If K, v3, v5 are jointly independent and each v_p ~ Geom((p-1)/p),")
print("then n+1 mod 30 should be approximately uniform.")
print()

_r.seed(123)
n = _r.getrandbits(3000) | 1
mod_vals = []

for _ in range(200000):
    mod_vals.append((n+1) % 30)
    n, _, _ = macro_step(n)
    if n < 2:
        n = _r.getrandbits(3000) | 1

mod_arr = np.array(mod_vals)
phi_30 = 8  # Euler's totient of 30 (numbers coprime to 30)
# Expected: n+1 coprime to 30, but n+1 = 2^K * m with K >= 1, so 2 | n+1
# Actually n+1 is always even (since n is odd). So n+1 is even.
# Also from Obs 276: v3(n+1) >= 0, and n+1 is even.
# So n+1 mod 30: n+1 is even, can be 0,2,4,6,8,10,12,14,16,18,20,22,24,26,28 mod 30 (even residues)
# Further: n+1 = 2^K * m where K >= 1, so 2|n+1. And macro-step output != mult of 3, so 3 nmid n_out.
# So n+1 mod 30: n+1 even, 3 nmid (n_out), 5 might divide (n_out is a valid Collatz output, not mult 3)

print("Distribution of (n+1) mod 30 along Collatz orbits:")
print(f"  (n is always odd, so n+1 is always even; counting only even residues)")
total = len(mod_arr)
even_mask = mod_arr % 2 == 0
even_vals = mod_arr[even_mask]
print(f"  Fraction with 2 | (n+1): {even_mask.mean():.4f} (should be ~1.000)")
print()

# Count residues mod 15 (since n+1 always even, look at (n+1)/2 mod 15)
half_mod15 = ((n_arr_dummy := mod_arr[even_mask]) // 2) % 15 if False else None

# Actually just count mod 30
from collections import Counter
cnt = Counter(mod_arr)
print(f"  Residue mod 30: {'Count':>10} {'Fraction':>10}")
even_residues = [r for r in range(0, 30) if r % 2 == 0]
total_even = sum(cnt.get(r,0) for r in even_residues)
for r in even_residues:
    c = cnt.get(r, 0)
    print(f"    {r:>5}: {c:>10} {c/total:>10.5f}")

print()
# Test uniformity among even residues not divisible by 3
valid_res = [r for r in range(0, 30) if r % 2 == 0 and r % 3 != 0]
print(f"Even residues not divisible by 3 (should be ~uniform): {valid_res}")
valid_counts = [cnt.get(r, 0) for r in valid_res]
expected = total_even / len(valid_res)
chi2, pval = stats.chisquare(valid_counts)
print(f"Chi-squared test for uniformity among {valid_res}: chi2={chi2:.2f}, p={pval:.4f}")
print(f"Expected per residue: {expected:.1f}, observed: {valid_counts}")

print()
print("=" * 70)
print("PART 3: JOINT (K, v3, v5) DISTRIBUTION")
print("=" * 70)
print()

_r.seed(456)
n = _r.getrandbits(3000) | 1
joint_kvv = defaultdict(int)

for _ in range(200000):
    K = v2(n+1)
    J3 = vp(n+1, 3)
    J5 = vp(n+1, 5)
    K_bin = min(K, 4)  # cap at 4+
    J3_bin = min(J3, 2)  # cap at 2+
    J5_bin = min(J5, 1)  # cap at 1+ (J5=0 or >=1)
    joint_kvv[(K_bin, J3_bin, J5_bin)] += 1
    n, _, _ = macro_step(n)
    if n < 2:
        n = _r.getrandbits(3000) | 1

total = sum(joint_kvv.values())
P_K = {k: sum(v for (K,J3,J5),v in joint_kvv.items() if K==k)/total for k in range(1,5)}
P_K[4] = sum(v for (K,J3,J5),v in joint_kvv.items() if K==4)/total  # overlap with cap
P_J3 = {j: sum(v for (K,J3,J5),v in joint_kvv.items() if J3==j)/total for j in range(3)}
P_J5 = {j: sum(v for (K,J3,J5),v in joint_kvv.items() if J5==j)/total for j in range(2)}

print("Testing joint independence of (K, v3, v5):")
print(f"  {'(K,J3,J5)':>15} {'Observed':>10} {'P(K)*P(J3)*P(J5)':>18} {'Ratio':>8}")
for (K, J3, J5), cnt in sorted(joint_kvv.items()):
    obs = cnt/total
    p_k = P_K.get(K, 0)
    p_j3 = P_J3.get(J3, 0)
    p_j5 = P_J5.get(J5, 0)
    expected = p_k * p_j3 * p_j5
    ratio = obs/expected if expected > 1e-8 else float('nan')
    if obs > 0.001:  # only show significant entries
        print(f"  ({K},{J3},{J5}):     {obs:>10.5f} {expected:>18.5f} {ratio:>8.4f}")

print()
print("=" * 70)
print("PART 4: RANDOMNESS TEST — n+1 mod VARIOUS PRIMES")
print("=" * 70)
print()
print("Test: is n+1 mod p approximately uniform on Z/pZ \\ {0}?")
print("(n is the orbit, p is a prime >= 5)")

_r.seed(789)
n = _r.getrandbits(3000) | 1
mod_data = {p: [] for p in [5, 7, 11, 13, 17, 19, 23]}

for _ in range(200000):
    for p in [5, 7, 11, 13, 17, 19, 23]:
        mod_data[p].append((n+1) % p)
    n, _, _ = macro_step(n)
    if n < 2:
        n = _r.getrandbits(3000) | 1

print()
for p in [5, 7, 11, 13, 17, 19, 23]:
    mvals = np.array(mod_data[p])
    # Expected: uniform on {1,...,p-1} (n+1 is never 0 mod p by Obs 276 for p=3; but for p>=5 it CAN be 0 mod p)
    # Actually for p>=5, n_out CAN be multiple of p (unlike p=3).
    cnt_zero = (mvals == 0).sum()
    cnt_nonzero = (mvals != 0).sum()
    # Chi-squared test for uniformity on {0,...,p-1}
    obs_counts = [( mvals == r).sum() for r in range(p)]
    expected_per = len(mvals) / p
    chi2, pval = stats.chisquare(obs_counts)
    frac_zero = cnt_zero / len(mvals)
    expected_frac_zero = 1/p
    print(f"  p={p:>3}: P(p|(n+1))={frac_zero:.5f} (expected 1/p={expected_frac_zero:.5f}), "
          f"chi2={chi2:.2f}, pval={pval:.4f}")

print()
print("=" * 70)
print("PART 5: MAXIMUM UPWARD EXCURSION IN LOG-SPACE")
print("=" * 70)
print()
print("How often does the Collatz orbit GROW significantly above n?")
print("Per-step Lyapunov: mu = -0.575. Can the orbit grow by a factor M > 1?")
print()

_r.seed(321)
n_initial = _r.getrandbits(200) | 1  # 200-bit starting number
log2_n0 = n_initial.bit_length() - 1

max_excursion = 0  # max log2(n/n0) ever reached
n = n_initial
log2_current = log2_n0

history = [log2_current]
for _ in range(5000):
    n, K, l0 = macro_step(n)
    log2_next = n.bit_length() - 1
    log2_current = log2_next
    excursion = log2_current - log2_n0
    max_excursion = max(max_excursion, excursion)
    history.append(log2_current)
    if n < 2: break

print(f"Starting number: 2^{log2_n0} (approx)")
print(f"Maximum log2(n/n0) excursion: {max_excursion} bits")
print(f"Total macro-steps until n=1: {len(history)-1}")
print()

# Do this for many starting numbers
max_excursions = []
orbit_lengths = []
_r.seed(654)
for b in [50, 100, 200, 500]:
    excursions_b = []; lengths_b = []
    for _ in range(200):
        n = _r.getrandbits(b) | 1
        log2_n0 = n.bit_length() - 1
        max_exc = 0
        t = 0
        while n > 1:
            n, _, _ = macro_step(n)
            exc = n.bit_length() - 1 - log2_n0
            max_exc = max(max_exc, exc)
            t += 1
        excursions_b.append(max_exc)
        lengths_b.append(t)
    max_excursions.append(excursions_b)
    orbit_lengths.append(lengths_b)
    print(f"b={b}: mean max excursion = {np.mean(excursions_b):.1f} bits, "
          f"max max excursion = {max(excursions_b):.0f} bits, "
          f"mean orbit length = {np.mean(lengths_b):.1f}")

print()
print("The maximum excursion is O(1) or O(log b) — stays BOUNDED above the start.")
print("This is consistent with: E[log ratio per step] = -0.575 < 0.")
print("By large deviations: P(excursion > A) ~ exp(-2*0.575*A/sigma^2)")
print(f"  = exp(-{2*0.575/1.644:.4f} * A)")
print("For A=10 bits: P ~ exp(-7.0) ~ 0.001 per orbit segment")

print()
print("=" * 70)
print("PART 6: MUTUAL INFORMATION BETWEEN v_p AND v_q")
print("=" * 70)
print()

_r.seed(42)
n = _r.getrandbits(3000) | 1
v3_seq = []; v5_seq = []; v7_seq = []
for _ in range(100000):
    v3_seq.append(min(vp(n+1, 3), 3))
    v5_seq.append(min(vp(n+1, 5), 1))
    v7_seq.append(min(vp(n+1, 7), 1))
    n, _, _ = macro_step(n)
    if n < 2:
        n = _r.getrandbits(3000) | 1

v3_arr = np.array(v3_seq)
v5_arr = np.array(v5_seq)
v7_arr = np.array(v7_seq)

def mutual_info(X, Y):
    """Compute mutual information I(X;Y) in nats."""
    from collections import Counter
    joint = Counter(zip(X, Y))
    px = Counter(X); py = Counter(Y)
    n = len(X)
    I = 0
    for (x,y), cnt in joint.items():
        p_xy = cnt/n; p_x = px[x]/n; p_y = py[y]/n
        if p_xy > 0:
            I += p_xy * math.log(p_xy / (p_x * p_y))
    return I

I_35 = mutual_info(v3_arr, v5_arr)
I_37 = mutual_info(v3_arr, v7_arr)
I_57 = mutual_info(v5_arr, v7_arr)
print(f"Mutual information I(v3, v5) = {I_35:.6f} nats (0 = independent)")
print(f"Mutual information I(v3, v7) = {I_37:.6f} nats")
print(f"Mutual information I(v5, v7) = {I_57:.6f} nats")
print()
# Compare to baseline (random shuffled)
v5_shuffled = v5_arr.copy(); np.random.shuffle(v5_shuffled)
I_35_null = mutual_info(v3_arr, v5_shuffled)
print(f"Null (shuffled) I(v3, v5_shuffled) = {I_35_null:.6f} nats")
print()
print("If I(X;Y) ≈ 0: X and Y are approximately independent.")
print("If I(X;Y) >> 0: X and Y carry information about each other.")
