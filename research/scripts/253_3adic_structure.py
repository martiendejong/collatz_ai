"""
253_3adic_structure.py
======================
Test the 3-adic structure hypothesis: does v2[s] differ systematically by s mod 3?

HYPOTHESIS: The doubling map s -> 2s+1 maps:
  s == 0 mod 3 -> 2s+1 == 1 mod 3  (SWAP)
  s == 1 mod 3 -> 2s+1 == 0 mod 3  (SWAP)
  s == 2 mod 3 -> 2s+1 == 2 mod 3  (FIXED)

If v2[s] has a systematic "level" based on s mod 3:
  Mean(v2 | s==0 mod 3) ≠ Mean(v2 | s==1 mod 3)

Then the 0↔1 swap creates anti-correlation:
  v2[s] large (s==0) => v2[2s+1] small (2s+1==1)
  v2[s] small (s==1) => v2[2s+1] large (2s+1==0)

Test 1: Is Mean(v2[s==0]) > Mean(v2[s==1]) or vice versa?
Test 2: Within each s mod 3 class, does the anti-correlation hold separately?
Test 3: What is the v2 mean pattern by s mod 3, and how does it propagate?

Also: check s mod 9, s mod 27 structure (multi-level 3-adic pattern).
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

def analyze_3adic(k, lam):
    v, Nl = run_kl(k, lam)
    v2 = v[2::3]  # v2[s] = v at (s, r=2), length Nl
    s = np.arange(Nl, dtype=np.int64)

    print(f"\nk={k}, lam={lam}, Nl={Nl}:")

    # === TEST 1: Mean v2 by s mod 3 ===
    for r in range(3):
        mask = (s % 3 == r)
        print(f"  Mean v2[s=={r} mod 3] = {float(np.mean(v2[mask])):.6f}  "
              f"Std = {float(np.std(v2[mask])):.6f}  n={mask.sum()}")

    # Ratio between classes
    m0 = float(np.mean(v2[s%3==0]))
    m1 = float(np.mean(v2[s%3==1]))
    m2 = float(np.mean(v2[s%3==2]))
    print(f"  Ratios: 0/1={m0/m1:.4f}  0/2={m0/m2:.4f}  1/2={m1/m2:.4f}")

    # === TEST 2: Anti-correlation by mod 3 class ===
    max_s = (Nl - 9) // 9
    s0 = np.arange(max_s, dtype=np.int64)
    v2_s = v2[s0]
    v2_2sp1 = v2[(2*s0 + 1) % Nl]

    # Split by s0 mod 3
    print(f"\n  Anti-correlation Corr(v2[2s+1], v2[s]) by s mod 3:")
    for r in range(3):
        mask = (s0 % 3 == r)
        if mask.sum() > 2:
            c = float(np.corrcoef(v2_s[mask], v2_2sp1[mask])[0,1])
            r2sp1 = (2*r+1) % 3  # which class 2s+1 falls in
            print(f"    s=={r} mod 3 (2s+1=={r2sp1} mod 3): Corr = {c:+.4f}  n={mask.sum()}")

    # Overall corr
    c_all = float(np.corrcoef(v2_s, v2_2sp1)[0,1])
    print(f"  Overall Corr(v2[2s+1], v2[s]) = {c_all:+.4f}")

    # === TEST 3: CONDITIONAL MEAN ===
    # Given s==0 mod 3 (expected large v2), is v2[2s+1] (with 2s+1==1) small?
    # Given s==1 mod 3 (expected small v2), is v2[2s+1] (with 2s+1==0) large?
    s0_mod0 = s0[s0 % 3 == 0]
    s0_mod1 = s0[s0 % 3 == 1]

    if len(s0_mod0) > 0 and len(s0_mod1) > 0:
        v2_s_mod0 = v2[s0_mod0]; v2_2sp1_mod0 = v2[(2*s0_mod0+1)%Nl]
        v2_s_mod1 = v2[s0_mod1]; v2_2sp1_mod1 = v2[(2*s0_mod1+1)%Nl]
        print(f"\n  s==0 mod 3: Mean v2[s]={float(np.mean(v2_s_mod0)):.4f}  "
              f"Mean v2[2s+1]={float(np.mean(v2_2sp1_mod0)):.4f}")
        print(f"  s==1 mod 3: Mean v2[s]={float(np.mean(v2_s_mod1)):.4f}  "
              f"Mean v2[2s+1]={float(np.mean(v2_2sp1_mod1)):.4f}")
        print(f"  SWAP: v2[2s+1] with 2s+1==1 < v2[2s+1] with 2s+1==0? "
              f"{float(np.mean(v2_2sp1_mod0)) < float(np.mean(v2_2sp1_mod1))}")

    # === TEST 4: MULTI-LEVEL 3-ADIC PATTERN ===
    # v2 mean by s mod 9
    print(f"\n  v2 mean by s mod 9 (multi-level pattern):")
    for r in range(9):
        mask = (s % 9 == r)
        if mask.sum() > 0:
            print(f"    s=={r:>2} mod 9: mean={float(np.mean(v2[mask])):.5f}  n={mask.sum()}")

    # === TEST 5: THEORETICAL PREDICTION ===
    # If v2[s] ≈ C * min(v2[(2s+1)//3], ...) and v2 = a[s mod 3] (three-level model),
    # Then: a[0] ≈ C * min(a[(2*0+1)//3], ...) = C * a[0]?... (2*0+1=1, 1//3=0, so a[0])
    # Wait: (2s+1)//3. For s==0: (2*0+1)//3 = 1//3 = 0. For s==1: (2*1+1)//3 = 3//3 = 1. For s==2: (2*2+1)//3 = 5//3 = 1.
    # Hmm, so v2[0] ≈ C * a[0], v2[1] ≈ C * a[1], v2[2] ≈ C * a[1] — doesn't create difference.

    # Let me instead look at v2[s] vs v2[(2s+1)%Nl] directly (the map s -> 2s+1)
    v2_map = v2[(2*s+1) % Nl]
    c_map = float(np.corrcoef(v2, v2_map)[0,1])
    print(f"\n  Corr(v2[s], v2[(2s+1)%Nl]) for ALL s = {c_map:+.4f}")
    print(f"  (map s -> 2s+1 over all s, not just subset)")

    # === TEST 6: The specific s-subset used in Script 251 ===
    # g_r2 = np.arange(Nl3)[np.arange(Nl3)%3==2], j_star = (4*g_r2)%Nl
    Nl3 = Nl // 3
    g_r2 = np.arange(Nl3)[np.arange(Nl3)%3==2]
    j_star = (4*g_r2) % Nl

    # Script 251's parent: v2[j_star]; child: v2[(j_star+Nl3)%Nl]
    # The 's_prime' = j_star//3 corresponds to "m0" in Script 251
    # j_star = 4g for g==2 mod 3. j_star == 2 mod 3 always.
    # The "child" 3s+2 structure: js1//3 = m0.
    # v2[j_star] is NOT v2[s_prime], it's v2 at the full j_star index.

    # So the actual s values used are j_star (for "parent"?) No...
    # In Script 251: v2_s1 = v2[js1] where js1 = (j_star+Nl3)%Nl.
    # And m0 = js1//3. The parent is NOT v2[j_star] but v2[m0].
    # The ACTUAL parent-child pair analyzed is (v2[3m+2], v2[m]) where m = m0 = js1//3.
    m0_arr = ((j_star + Nl3) % Nl) // 3
    v2_3m2 = v2[(j_star + Nl3) % Nl]  # = v2[3*m0+2]
    v2_m = v2[m0_arr % Nl]             # = v2[m0]

    corr_pc = float(np.corrcoef(v2_3m2, v2_m)[0,1])
    print(f"\n  Script 251 parent-child Corr(v2[3m+2], v2[m]) = {corr_pc:+.4f}")

    # m0 values: m0 = (4g//3 + Nl9) % (Nl//3) roughly
    # What mod 3 are 3m+2 and m?
    # 3m+2 == 2 mod 3 ALWAYS
    # m can be any value mod 3
    print(f"  m0 mod 3 distribution: {dict(zip(*np.unique(m0_arr%3, return_counts=True)))}")
    print(f"  3m+2 mod 3 = {dict(zip(*np.unique((3*m0_arr+2)%3, return_counts=True)))}")

    # The key: v2[3m+2] is always at position ==2 mod 3.
    # v2[m] is at position =={0,1,2} mod 3 depending on m.
    # If v2[s==2] and v2[s==0,1] have different means, this contributes to the correlation.
    print(f"\n  v2[m0] by m0 mod 3:")
    for r in range(3):
        mask = (m0_arr % 3 == r)
        if mask.sum() > 0:
            print(f"    m0=={r} mod 3: mean v2[m0]={float(np.mean(v2_m[mask])):.5f}  "
                  f"mean v2[3m+2]={float(np.mean(v2_3m2[mask])):.5f}  n={mask.sum()}")

    return c_map

print("253: 3-adic structure of K-L eigenvector")
print("Testing: does s mod 3 pattern explain Corr(v2[2s+1], v2[s]) < 0?")
print("="*70)

analyze_3adic(8, 1.70)

print("\n\n=== Lambda scan k=8: corr over ALL s (map s->2s+1) ===")
for lam in [1.30, 1.50, 1.70, 1.90, 2.00]:
    v, Nl = run_kl(8, lam)
    v2 = v[2::3]
    s = np.arange(Nl, dtype=np.int64)
    c = float(np.corrcoef(v2, v2[(2*s+1)%Nl])[0,1])
    means = [float(np.mean(v2[s%3==r])) for r in range(3)]
    print(f"  lam={lam:.2f}: Corr(v2, v2[2s+1])={c:+.4f}  "
          f"means[0]={means[0]:.4f} [1]={means[1]:.4f} [2]={means[2]:.4f}")
    sys.stdout.flush()

print("\ndone")
