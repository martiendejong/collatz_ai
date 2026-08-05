"""
255_recursive_anticorr.py
==========================
Verify that the anti-correlation Corr(v2[2s+1], v2[s]) has a RECURSIVE structure:

Level 1 (mod 3 between-class):
  a0 < a1 + map 0<->1 swap => anti-corr
  BUT: at large lambda, a0/a1 -> 1, between-class effect vanishes

Level 2 (mod 9 within class 2->2):
  The map s->2s+1 on the s==2 mod 3 class maps:
  s==2 mod 9 -> 2s+1==5 mod 9 (and 5->2, cycle of length 2)
  s==8 mod 9 -> 2s+1==8 mod 9 (fixed point)
  Means: 0.428, 0.198, 0.462 (for 2,5,8 mod 9)
  So within class 2: 2<->5 swap with mean(2)=0.428 >> mean(5)=0.198 => anti-corr

Level 2 (mod 9 within class 0->1 and 1->0):
  These are already between-class at the mod-3 level, but within each:
  s==0 mod 9: 2s+1==1 mod 9
  s==3 mod 9: 2s+1==7 mod 9
  s==6 mod 9: 2s+1==4 mod 9
  Means at mod 9: 0.167, 0.109, 0.144 (for s mod 9 = 0,3,6)
                  0.236, 0.432, 0.316 (for 2s+1 mod 9 = 1,7,4)
  Pairs: (0.167, 0.236), (0.109, 0.432), (0.144, 0.316)
  Each pair has the TARGET (2s+1 class) > SOURCE (s class)? Partially.

HYPOTHESIS: the anti-correlation is SELF-SIMILAR.
  At level m (mod 3^m), the mean structure follows the SAME pattern as at level 1,
  just with means averaged over sub-cells.

If true: Corr(v2[2s+1], v2[s]) can be decomposed as:
  Level 1 contribution: from mod-3 class means
  Level 2 contribution: within-class mod-9 sub-structure
  Level m contribution: from mod-3^m structure
  Total = sum of all levels

TESTS:
1. Compute Corr(v2[2s+1], v2[s]) residual after removing mod-3 mean
2. Compute whether residual has a mod-9 structure
3. Recursively: compute how many levels are needed to explain the full anti-corr
4. At large lambda (where between-class fails), show within-class is the driver
"""
import numpy as np
from math import log2
import sys

ALPHA = log2(3.0)
N_ITER = 600

def run_kl(k, lam, n_iter=N_ITER):
    A  = lam ** -2.0
    B1 = lam ** (ALPHA - 2.0)
    B3 = lam ** (ALPHA - 1.0)
    N  = 3 ** (k - 1)
    Nl = N // 3
    i  = np.arange(N, dtype=np.int64)
    T4 = (4 * i + 2) % N
    s_arr, r_arr = np.divmod(i, 3)
    m0, m2 = (r_arr == 0), (r_arr == 2)
    R1 = (4 * s_arr) % Nl
    R3 = (2 * s_arr + 1) % Nl
    v = np.ones(N, dtype=np.float64)
    for _ in range(n_iter):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        w  = A * v[T4]
        w[m2] += B3 * cb[R3[m2]]
        w[m0] += B1 * cb[R1[m0]]
        v = w / w.max()
    return v, Nl

def multilevel_corr_decomp(k, lam):
    """Decompose Corr(v2[2s+1], v2[s]) by progressively including more 3-adic levels."""
    v, Nl = run_kl(k, lam)
    v2 = v[2::3]
    s = np.arange(Nl, dtype=np.int64)

    # Overall correlation
    max_s = (Nl - 9) // 9
    s0 = np.arange(max_s, dtype=np.int64)
    v2_s = v2[s0]
    v2_2sp1 = v2[(2*s0+1) % Nl]
    corr_actual = float(np.corrcoef(v2_s, v2_2sp1)[0,1])

    print(f"\nk={k}, lam={lam:.2f}: Actual Corr(v2[2s+1], v2[s]) = {corr_actual:+.4f}")

    results = []
    # Level m: replace v2[s] by its mod-3^m class mean, compute implied correlation
    for m in range(1, 7):
        mod_m = 3 ** m
        # Class means
        classes_s = s0 % mod_m
        classes_2sp1 = (2*s0+1) % mod_m
        unique_s = np.unique(classes_s)

        # Compute mean of v2 for each class (using all Nl positions, not just s0)
        # because s0 is only (Nl-9)//9 elements
        all_s = np.arange(Nl, dtype=np.int64)
        mean_map = np.zeros(mod_m, dtype=np.float64)
        for c in range(mod_m):
            mask = (all_s % mod_m == c)
            if mask.sum() > 0:
                mean_map[c] = float(np.mean(v2[mask]))

        # Approximate v2 by class means
        v2_proxy = mean_map[s0 % mod_m]
        v2_2sp1_proxy = mean_map[(2*s0+1) % mod_m]

        # Correlation of class-mean proxies
        if np.std(v2_proxy) > 0 and np.std(v2_2sp1_proxy) > 0:
            c_proxy = float(np.corrcoef(v2_proxy, v2_2sp1_proxy)[0,1])
        else:
            c_proxy = 0.0

        results.append((mod_m, c_proxy))
        print(f"  Level m={m} (mod {mod_m:>5}): proxy Corr = {c_proxy:+.4f}  "
              f"fraction of total = {c_proxy/corr_actual:+.4f}")

    # Now: partition the actual corr into levels
    # Level 1 contribution = between-class at mod 3
    # Level 2 contribution = between-class at mod 9 MINUS level 1 (additional explained)
    print(f"\n  Incremental level contributions:")
    prev = 0.0
    for mod_m, c_proxy in results:
        increment = c_proxy - prev
        print(f"    mod {mod_m:>5}: additional Corr = {increment:+.4f}  cumulative = {c_proxy:+.4f}")
        prev = c_proxy
    print(f"    Unexplained (within mod {results[-1][0]}): {corr_actual - prev:+.4f}")

    return corr_actual, results

def within_class_2_recursion(k, lam):
    """Verify that within class 2->2, the anti-corr has a mod-9 sub-structure."""
    v, Nl = run_kl(k, lam)
    v2 = v[2::3]
    s = np.arange(Nl, dtype=np.int64)

    max_s = (Nl - 9) // 9
    s0 = np.arange(max_s, dtype=np.int64)

    # Extract only s==2 mod 3
    mask_2 = (s0 % 3 == 2)
    s2 = s0[mask_2]  # s values with s==2 mod 3
    # 2s+1 when s==2 mod 3: 2*2+1=5==2 mod 3, confirmed stays in class 2
    v2_s2 = v2[s2]
    v2_2sp1_s2 = v2[(2*s2+1) % Nl]
    corr_22 = float(np.corrcoef(v2_s2, v2_2sp1_s2)[0,1])
    print(f"\n  Within class 2->2 (s==2 mod 3): Corr = {corr_22:+.4f}  n={len(s2)}")

    # What is 2s+1 mod 9 when s==2 mod 3? (s can be 2,5,8 mod 9)
    print(f"  Mod-9 structure within class 2:")
    print(f"    s==2 mod 9 -> 2s+1==5 mod 9")
    print(f"    s==5 mod 9 -> 2s+1==2 mod 9  (swap with 2)")
    print(f"    s==8 mod 9 -> 2s+1==8 mod 9  (fixed)")

    all_s = np.arange(Nl, dtype=np.int64)
    print(f"  Means by class:")
    for r in range(9):
        mask = (all_s % 9 == r)
        if mask.sum() > 0:
            print(f"    s=={r} mod 9: mean = {float(np.mean(v2[mask])):.5f}")

    # Between-class contribution within class 2 (using mod-9 means)
    # The map 2->5->2 creates anti-corr if mean(class 2 mod 9) != mean(class 5 mod 9)
    mean9 = np.zeros(9)
    for r in range(9):
        mask = (all_s % 9 == r)
        if mask.sum() > 0:
            mean9[r] = float(np.mean(v2[mask]))

    # Pairs in class 2 (mod 3): s in {2,5,8} mod 9
    # Map: 2->5, 5->2, 8->8
    # Only consider s2 with valid mod-9 sub-class
    # Compute proxy correlation using mod-9 means within class 2
    s2_9 = s2 % 9  # s mod 9 for s in class 2
    proxy_v2 = mean9[s2_9]
    proxy_v2_2sp1 = mean9[(2*s2+1) % 9]
    if np.std(proxy_v2) > 0 and np.std(proxy_v2_2sp1) > 0:
        c_proxy_22 = float(np.corrcoef(proxy_v2, proxy_v2_2sp1)[0,1])
        print(f"\n  Within class 2, between-subclass (mod-9) Corr = {c_proxy_22:+.4f}")
        print(f"  Fraction of Corr_22 explained by mod-9: {c_proxy_22/corr_22:.4f}")
        # Means by s mod 9 within class 2
        for r9 in [2, 5, 8]:
            mask9 = (s2 % 9 == r9)
            n = mask9.sum()
            r9_2sp1 = (2*r9+1) % 9
            print(f"    s=={r9} mod 9 -> 2s+1=={r9_2sp1} mod 9: "
                  f"mean(v2[s])={mean9[r9]:.4f} mean(v2[2s+1])={mean9[r9_2sp1]:.4f}  n={n}")

print("255: Recursive multi-level decomposition of anti-correlation")
print("="*70)

k, lam = 8, 1.70
multilevel_corr_decomp(k, lam)
within_class_2_recursion(k, lam)

print("\n\n=== Lambda scan: how many levels needed? ===")
for lam in [1.30, 1.50, 1.70, 1.90, 2.00]:
    v, Nl = run_kl(8, lam)
    v2 = v[2::3]
    s = np.arange(Nl, dtype=np.int64)
    max_s = (Nl - 9) // 9
    s0 = np.arange(max_s, dtype=np.int64)
    corr_actual = float(np.corrcoef(v2[s0], v2[(2*s0+1)%Nl])[0,1])

    parts = []
    for m in range(1, 7):
        mod_m = 3 ** m
        all_s = np.arange(Nl, dtype=np.int64)
        mean_map = np.array([float(np.mean(v2[all_s % mod_m == c])) for c in range(mod_m)])
        v2_proxy = mean_map[s0 % mod_m]
        v2_2sp1_proxy = mean_map[(2*s0+1) % mod_m]
        if np.std(v2_proxy) > 0 and np.std(v2_2sp1_proxy) > 0:
            c_proxy = float(np.corrcoef(v2_proxy, v2_2sp1_proxy)[0,1])
        else:
            c_proxy = 0.0
        parts.append(c_proxy)

    print(f"lam={lam:.2f}: actual={corr_actual:+.4f}  "
          f"L1={parts[0]:+.4f}  L2={parts[1]:+.4f}  L3={parts[2]:+.4f}  "
          f"L4={parts[3]:+.4f}  L5={parts[4]:+.4f}  L6={parts[5]:+.4f}")
    sys.stdout.flush()

print("\ndone")
