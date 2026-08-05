"""
274_pointwise_cov2_analysis.py
==============================
POINTWISE per-column CoV² comparison between v0 and v2 column triplets.

Goal: determine whether CoV²(v2-col) > CoV²(v0-col) holds POINTWISE
(per individual column) or only ON AVERAGE across sub-groups.

The Jensen monotonicity theorem gives m2m(v2-col) < m2m(v0-col) IF
CoV²(v2-col) > CoV²(v0-col) FOR EACH COLUMN (pointwise). If the
inequality only holds on average, Jensen does not directly apply.

Analysis plan:
  1. Per-column pointwise fraction: CoV²(v2-col[j3]) > CoV²(v0-col[j3]).
  2. Sub-group fractions (r=0,1,2).
  3. For r=2: ANALYTICAL argument via sub-class structure.
     Condition: lambda^2*(1+t^4) > 1+t^2 -- derive exactly.
  4. Conditional m2m check: columns where CoV²(v2)>CoV²(v0): do they
     always have m2m(v2)<m2m(v0)?
  5. NET SIGNED EFFECT: is the weighted sum always positive?

NEW ANALYTICAL RESULT (derived from Obs 471 global equations):
  Var(v2-col) > Var(v0-col) iff lambda^2*(1+t^4) > 1+t^2.
  Proof: the condition simplifies to (1-t^2)(lambda^2 - 1 - t^2) > 0,
  which holds iff lambda^2 > 1+t^2.
  For lambda=2: always holds (4 > 2). For lambda >= 1.10: holds when
  t < sqrt(lambda^2 - 1).
"""
import numpy as np
from math import log2, sqrt
import sys

ALPHA = log2(3.0)

def run_kl(k, lam, n_iter=800):
    A  = lam**-2.0; B1 = lam**(ALPHA-2.0); B3 = lam**(ALPHA-1.0)
    N  = 3**(k-1); Nl = N//3
    i  = np.arange(N, dtype=np.int64)
    T4 = (4*i+2) % N
    s_arr, r_arr = np.divmod(i, 3)
    m0_mask, m2_mask = (r_arr==0), (r_arr==2)
    R1 = (4*s_arr) % Nl; R3 = (2*s_arr+1) % Nl
    v = np.ones(N, dtype=np.float64)
    for _ in range(n_iter):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        w  = A*v[T4]
        w[m2_mask] += B3*cb[R3[m2_mask]]
        w[m0_mask] += B1*cb[R1[m0_mask]]
        v = w/w.max()
    return v, Nl, A, B1, B3

def compute_R(t, lam):
    return (t**2 + lam) / (1 + t*lam)

def compute_t(v, Nl, A, rho):
    """Estimate t = A/rho from Perron eigenvalue."""
    return A / rho

def col_cov2(col_matrix):
    """col_matrix: (n_cols, 3). Returns per-column CoV²."""
    m = col_matrix.mean(axis=1, keepdims=True)
    v = ((col_matrix - m)**2).mean(axis=1)
    return v / (m.squeeze()**2 + 1e-300)

def col_m2m(col_matrix):
    """min-to-mean ratio per column."""
    return col_matrix.min(axis=1) / col_matrix.mean(axis=1)

print("="*70)
print("SCRIPT 274: POINTWISE CoV² COMPARISON v0 vs v2 COLUMN TRIPLETS")
print("="*70)

# ======================================================================
# SECTION 1: ANALYTICAL CONDITION
# ======================================================================
print()
print("SECTION 1: ANALYTICAL CONDITION lambda^2 > 1+t^2")
print()
print("From Obs 471 equations:")
print("  P = A'*Q + B'    (v0-col variance equation)")
print("  Q = t^2*A'*P + lam^2*B'  (v2-col variance equation)")
print("  where A' = (A/rho)^2 = t^2, B' = (B1/rho)^2 * C_cb")
print()
print("Var(v2-col) > Var(v0-col) iff:")
print("  (lam^2-1)*B' > t^2*(Q - t^2*P)")
print("  Using B' = P - t^2*Q and Q/P = (t^4+lam^2)/(1+t^2*lam^2):")
print("  Simplifies to: (1-t^2)*(lam^2 - 1 - t^2) > 0")
print("  => lambda^2 > 1+t^2")
print()
print("Maximum of 1+t^2 over t in (0,1): approaches 2 as t->1.")
print("But as lam->1+, t->1-. So at lam=1: need lam^2 > 2 (fails).")
print("Critical: for what (lam,t) pairs does lam^2 > 1+t^2 hold?")
print()

# For each tested lam, compute the critical t threshold
for lam in [1.10, 1.20, 1.30, 1.50, 1.70, 1.90, 2.00]:
    t_crit = sqrt(max(0, lam**2 - 1))
    print(f"  lam={lam:.2f}: lam^2={lam**2:.4f}, t_crit=sqrt(lam^2-1)={t_crit:.4f}")
    print(f"    => Var(v2)>Var(v0) when t < {t_crit:.4f}")

print()
print("Now check actual t values for our test cases:")
print()

# ======================================================================
# SECTION 2: COMPUTE ACTUAL t, R, AND VERIFY CONDITION
# ======================================================================
print("SECTION 2: ACTUAL (lam, t) VALUES AND POINTWISE FRACTION")
print()

test_cases = [
    (1.30, 8), (1.30, 10), (1.30, 12),
    (1.50, 8), (1.50, 10), (1.50, 12),
    (1.70, 8), (1.70, 10), (1.70, 12),
    (1.90, 8), (1.90, 10), (1.90, 12),
    (2.00, 8), (2.00, 10),
]

results = []

for lam, k in test_cases:
    v, Nl, A, B1, B3 = run_kl(k, lam)
    Nl3 = Nl // 3
    j3 = np.arange(Nl3, dtype=np.int64)

    v0 = v[0::3]; v1 = v[1::3]; v2 = v[2::3]

    # Estimate rho from v0/v2 ratio (eigenvalue property: v_next = rho * v_current scaled)
    # Better: t = A/rho where rho satisfies v0.mean() / v2.mean() = rho / (A*t + B1*c0/v0.mean())
    # Simplest: use the R formula. R = mean_v2/mean_v0 = (t^2+lam)/(1+t*lam)
    # From this: t^2 - R*t*lam + lam - R = 0... solve for t.
    R_actual = float(v2.mean() / v0.mean())
    # Solve: (t^2+lam)/(1+t*lam) = R => t^2 + lam = R + R*t*lam
    # t^2 - R*lam*t + (lam-R) = 0
    disc = (R_actual*lam)**2 - 4*(lam - R_actual)
    if disc >= 0:
        t1 = (R_actual*lam - sqrt(disc)) / 2
        t2 = (R_actual*lam + sqrt(disc)) / 2
        # t should be in (0,1)
        t_val = t1 if 0 < t1 < 1 else t2
    else:
        t_val = float('nan')

    cond_holds = lam**2 > 1 + t_val**2
    t_crit = sqrt(max(0, lam**2 - 1))

    # Per-column CoV² and m2m
    col_v0 = np.stack([v0[j3], v0[j3+Nl3], v0[j3+2*Nl3]], axis=1)
    col_v2 = np.stack([v2[j3], v2[j3+Nl3], v2[j3+2*Nl3]], axis=1)

    cov2_v0 = col_cov2(col_v0)
    cov2_v2 = col_cov2(col_v2)
    m2m_v0  = col_m2m(col_v0)
    m2m_v2  = col_m2m(col_v2)

    frac_cov2_correct = float(np.mean(cov2_v2 > cov2_v0))
    frac_m2m_correct  = float(np.mean(m2m_v2 < m2m_v0))
    global_m2m_gap    = float(np.mean(m2m_v0 - m2m_v2))

    results.append({
        'lam': lam, 'k': k, 't': t_val, 'R': R_actual,
        'cond_holds': cond_holds, 't_crit': t_crit,
        'frac_cov2': frac_cov2_correct, 'frac_m2m': frac_m2m_correct,
        'global_gap': global_m2m_gap,
        'col_v0': col_v0, 'col_v2': col_v2,
        'cov2_v0': cov2_v0, 'cov2_v2': cov2_v2,
        'm2m_v0': m2m_v0, 'm2m_v2': m2m_v2,
        'j3': j3
    })

    print(f"lam={lam:.2f} k={k:2d}: t={t_val:.4f} t_crit={t_crit:.4f} "
          f"cond(lam^2>1+t^2)={cond_holds} "
          f"frac(CoV2_v2>CoV2_v0)={frac_cov2_correct:.3f} "
          f"frac(m2m_v2<m2m_v0)={frac_m2m_correct:.3f} "
          f"global_gap={global_m2m_gap:.5f}")

print()

# ======================================================================
# SECTION 3: SUB-GROUP POINTWISE ANALYSIS
# ======================================================================
print("="*70)
print("SECTION 3: SUB-GROUP POINTWISE FRACTIONS (lam=1.70, k=10)")
print()

r70 = [r for r in results if r['lam'] == 1.70 and r['k'] == 10][0]
j3 = r70['j3']
cov2_v0 = r70['cov2_v0']
cov2_v2 = r70['cov2_v2']
m2m_v0  = r70['m2m_v0']
m2m_v2  = r70['m2m_v2']
Nl3 = len(j3)

for r_grp in [0, 1, 2]:
    mask = (j3 % 3 == r_grp)
    n = mask.sum()
    cv0 = cov2_v0[mask]; cv2 = cov2_v2[mask]
    mv0 = m2m_v0[mask];  mv2 = m2m_v2[mask]

    frac_c = float(np.mean(cv2 > cv0))
    frac_m = float(np.mean(mv2 < mv0))
    avg_cov2_diff = float(np.mean(cv2 - cv0))
    avg_m2m_diff  = float(np.mean(mv0 - mv2))

    print(f"r={r_grp} (n={n}):")
    print(f"  frac(CoV2_v2>CoV2_v0) = {frac_c:.4f}")
    print(f"  frac(m2m_v2<m2m_v0)  = {frac_m:.4f}")
    print(f"  mean(CoV2_v2-CoV2_v0) = {avg_cov2_diff:.6f}")
    print(f"  mean(m2m_v0-m2m_v2)  = {avg_m2m_diff:.6f}")

    # Correlation: does higher CoV2 diff => higher m2m_v0-m2m_v2?
    cov2_diff = cv2 - cv0
    m2m_diff  = mv0 - mv2
    if cov2_diff.std() > 1e-10:
        corr = float(np.corrcoef(cov2_diff, m2m_diff)[0,1])
        print(f"  corr(CoV2_diff, m2m_diff) = {corr:.4f}")
    print()

# ======================================================================
# SECTION 4: CONDITIONAL CHECK - when CoV2_v2 > CoV2_v0, does m2m_v2 < m2m_v0?
# ======================================================================
print("="*70)
print("SECTION 4: CONDITIONAL CHECK across all test cases")
print()
print("Question: when CoV2(v2-col) > CoV2(v0-col) pointwise,")
print("          is m2m(v2-col) < m2m(v0-col) for that column?")
print()

for res in results[:6]:  # First 6 test cases
    lam, k = res['lam'], res['k']
    cov2_diff = res['cov2_v2'] - res['cov2_v0']
    m2m_diff  = res['m2m_v0'] - res['m2m_v2']
    # When CoV2 diff > 0, is m2m diff > 0?
    pos_cov2_mask = cov2_diff > 0
    neg_cov2_mask = cov2_diff <= 0
    if pos_cov2_mask.sum() > 0:
        frac_given_pos = float(np.mean(m2m_diff[pos_cov2_mask] > 0))
        avg_m2m_given_pos = float(np.mean(m2m_diff[pos_cov2_mask]))
    else:
        frac_given_pos = float('nan')
        avg_m2m_given_pos = float('nan')
    if neg_cov2_mask.sum() > 0:
        frac_given_neg = float(np.mean(m2m_diff[neg_cov2_mask] > 0))
        avg_m2m_given_neg = float(np.mean(m2m_diff[neg_cov2_mask]))
    else:
        frac_given_neg = float('nan')
        avg_m2m_given_neg = float('nan')
    n_pos = pos_cov2_mask.sum()
    n_neg = neg_cov2_mask.sum()
    print(f"lam={lam:.2f} k={k:2d}:")
    print(f"  Given CoV2_v2>CoV2_v0 (n={n_pos}): "
          f"frac(m2m_v2<m2m_v0)={frac_given_pos:.3f}, "
          f"avg_m2m_gap={avg_m2m_given_pos:.5f}")
    print(f"  Given CoV2_v2<=CoV2_v0 (n={n_neg}): "
          f"frac(m2m_v2<m2m_v0)={frac_given_neg:.3f}, "
          f"avg_m2m_gap={avg_m2m_given_neg:.5f}")

print()

# ======================================================================
# SECTION 5: WEIGHTED NET EFFECT ANALYSIS
# ======================================================================
print("="*70)
print("SECTION 5: NET SIGNED EFFECT — can we bound the aggregate?")
print()
print("If we CANNOT prove pointwise CoV2_v2 > CoV2_v0,")
print("can we at least bound the net m2m difference analytically?")
print()
print("Net effect = sum_j [m2m_v0(j) - m2m_v2(j)]")
print("           = m2m_v0(global) - m2m_v2(global) * (n_cols)")
print("           = E[m2m_v0] - E[m2m_v2] > 0 iff step (3b) holds.")
print()

# Compute by sub-group weighted sum
lam, k = 1.70, 10
r70 = [r for r in results if r['lam'] == lam and r['k'] == k][0]
j3 = r70['j3']
m2m_v0 = r70['m2m_v0']; m2m_v2 = r70['m2m_v2']

total_gap = float(np.mean(m2m_v0 - m2m_v2))
print(f"lam={lam} k={k}: total gap E[m2m_v0-m2m_v2] = {total_gap:.6f}")
print()

for r_grp in [0, 1, 2]:
    mask = (j3 % 3 == r_grp)
    n = mask.sum()
    gap_r = float(np.mean((m2m_v0 - m2m_v2)[mask]))
    contribution = float(n * gap_r / len(j3))
    print(f"  r={r_grp}: mean_gap={gap_r:.6f}, n={n}, contribution={contribution:.6f}")

print(f"  Sum of contributions = {sum([float(np.mean((m2m_v0 - m2m_v2)[j3%3==r])) / 3 for r in [0,1,2]]):.6f}")
print()

# ======================================================================
# SECTION 6: ANALYTICAL CLOSURE ATTEMPT FOR r=2
# ======================================================================
print("="*70)
print("SECTION 6: ANALYTICAL CLOSURE FOR r=2 SUB-GROUP")
print()
print("For r=2: BOTH v0 and v2 use sc2-cb (same cb class).")
print("The ONLY difference is the Z-input source:")
print("  v0-col r=2: Z from {v2[sc2]}")
print("  v2-col r=2: Z from {t * v0[sc1]}")
print()
print("Within-col variance of Z-inputs:")
print("  v0: var({v2[T4(j3+s*Nl3)]}) -- all in sc2 of v2")
print("  v2: t^2 * var({v0[sigma1(T4(j3+s*Nl3))]}) -- all in sc1 of v0")
print()
print("If the WITHIN-COL variance of sc2-v2 > t^2 * within-col var of sc1-v0,")
print("then CoV2(v2-col r=2) > CoV2(v0-col r=2) for EACH column in r=2.")
print()

# For r=2 sub-group, compute per-column within-col variance of Z-inputs
lam, k = 1.70, 10
v, Nl, A, B1, B3 = run_kl(k, lam)
Nl3 = Nl // 3
j3 = np.arange(Nl3, dtype=np.int64)
v0 = v[0::3]; v1 = v[1::3]; v2 = v[2::3]

R_actual = float(v2.mean() / v0.mean())
disc = (R_actual*lam)**2 - 4*(lam - R_actual)
t_val = (R_actual*lam - sqrt(disc)) / 2
rho_est = A / t_val

rho_A = A / rho_est; rho_B1 = B1 / rho_est; rho_B3 = B3 / rho_est

jj_r2 = j3[j3 % 3 == 2]
print(f"lam={lam}, k={k}: t={t_val:.4f}, R={R_actual:.5f}")
print(f"  n(r=2) = {len(jj_r2)}")
print()

# v0-col r=2: Z-inputs are v2[T4(j3+s*Nl3)] = v2[4*(j3+s*Nl3) % Nl]
# For j3 = 3*m+2 (3-indexed): T4(j3+s*Nl3) = (4*(j3+s*Nl3)+2) % N in full array
# but in Nl-indexed: the v2 sub-index of T4 is ((4*(j3+s*Nl3)+2)-2)/3 = 4*(j3+s*Nl3)/3...
# Actually let me compute directly from the eigenvector.

# For each r=2 column, the three Z-inputs to v0-col:
# v0[j] = A/rho * v2[T4(j)] + B1/rho * cb[R1(j)]  (for r=1 mode, i.e., actual j%3=0 in v0 indexing)
# Wait -- need to re-read the operator structure.

# K-L operator: v_next[i] = A*v[T4(i)] + Br*cb[Rr(i)] for r = i%3.
# For the eigenvector: v0 = v[0::3], and in FULL ARRAY index, v0-elements are {3*j3, j3 in 0..Nl-1}.
# For element at full-index i=3*j3 (r=0 in full): r=0 mod 3 = 0.
# The cb contribution uses B1 (r=0 uses R1).
# But WAIT: column-triplets group v0[j3], v0[j3+Nl3], v0[j3+2*Nl3] in the v0-sub-array.
# In the FULL array, these are at indices 3*j3, 3*(j3+Nl3), 3*(j3+2*Nl3).
# All have full_r = 0 (they're in the v0 sub-array).
# So the K-L operator for ALL v0 elements uses R1, B1.

# The v0-sub-array j3 has r_sub = j3 % 3. This is the sub-class within v0.
# For r_sub=2 columns of v0:
# Z-input: v2[T4(3*(j3+s*Nl3)) // 3] = v2[(4*3*(j3+s*Nl3)+2) // 3 % Nl]

# Exact computation:
# T4_full(i) = (4*i + 2) % N where N = 3^(k-1)
# For i = 3*(j3+s*Nl3):
#   T4_full = (4*3*(j3+s*Nl3) + 2) % N = (12*(j3+s*Nl3) + 2) % N
# The v2-sub-index of this: since T4 gives an r=2 full element (as proved in Script 273):
#   v2_subindex = (T4_full - 2) / 3 = (12*(j3+s*Nl3)+2-2)/3 = 4*(j3+s*Nl3)
#   Reduced mod Nl: 4*(j3+s*Nl3) % Nl = 4*j3 % Nl (since 4*Nl3*s % Nl = 0 when Nl3=Nl//3)
# Hmm: 4*Nl3 % Nl = 4*(Nl//3) % Nl. For Nl = 3^(k-2): Nl3 = 3^(k-3).
#   4*Nl3 % Nl = 4*3^(k-3) % 3^(k-2) = 4*3^(k-3) % 3^(k-2).
#   If k>=3: 3^(k-2) = 3*3^(k-3). So 4*3^(k-3) mod 3*3^(k-3) = 3^(k-3).
#   So the three Z-inputs for s=0,1,2 are at v2-indices:
#     4*j3 % Nl, (4*j3 + Nl3) % Nl, (4*j3 + 2*Nl3) % Nl
# These are THREE DIFFERENT ELEMENTS of v2 separated by Nl3!
# Their sc-types: (4*j3) % 3, (4*j3+Nl3) % 3, (4*j3+2*Nl3) % 3.
# Since Nl3 % 3 = 0: all three have sc-type = (4*j3) % 3 = j3 % 3 = 2.
# So yes, all three Z-inputs are from sc2 of v2.
# But crucially: they are at DIFFERENT POSITIONS within sc2!
# Specifically: 4*j3, 4*j3+Nl3, 4*j3+2*Nl3 (mod Nl) within the Nl3-indexed sub-class.

# So the WITHIN-COL variance of Z-inputs for v0-col r=2 is:
# var(v2[4*j3 % Nl], v2[(4*j3+Nl3) % Nl], v2[(4*j3+2*Nl3) % Nl])
# These are elements at positions {4*j3, 4*j3+Nl3, 4*j3+2*Nl3} within v2.
# In sc2-of-v2 sub-index: {4*j3, 4*j3+Nl3, 4*j3+2*Nl3} all mod Nl.
# Equivalently, in the COLUMN TRIPLET: these are {v2[alpha], v2[alpha+Nl3], v2[alpha+2*Nl3]}
# where alpha = (4*j3) % Nl. This is ANOTHER COLUMN TRIPLET of v2 (with different j3)!

# KEY INSIGHT: The Z-inputs for v0-col-r2[j3] form a COLUMN-TRIPLET of v2 at j3' = (4*j3 % Nl3).
# And the v2 column-triplet at j3' is ALSO a column with j3'%3 = (4*j3)%3 % 3...

# This establishes a MAP between v0-col-r2 columns and v2-col columns.
# The within-col variance of the v0-col-r2 Z-inputs = variance of the v2-col at j3' = 4*j3 % Nl3.

# Similarly for v2-col r=2:
# Z-inputs are t * v0[sigma1(T4(3*(j3+2)+s*Nl3))]...
# Let me compute directly.

# Direct numerical approach: compute per-column within-col Z-variance
# for both v0-col-r2 and v2-col-r2.

# For v0-col-r2: Z-inputs are v2 at positions 4*j3, 4*j3+Nl3, 4*j3+2*Nl3 (mod Nl)
Z0_for_v0r2 = np.stack([
    v2[(4*jj_r2) % Nl],
    v2[(4*jj_r2 + Nl3) % Nl],
    v2[(4*jj_r2 + 2*Nl3) % Nl]
], axis=1)

# For v2-col-r2: v2[j] = A/rho * v1[T4(j)] + B3/rho * cb[R3(j)]
# v1 at sub-class 2: v1[j3_sub] = (A/rho)*v0[sigma1(j3_sub)] + 0 (only for r=1,2 types)
# Actually v1[sc2] uses v0 at sc1 of v0.
# Let's compute directly: for v2-col at jj_r2,
# The T4 mapping for v2 elements at full index 3*jj_r2 + 2 (since v2 = v[2::3]):
# Full index of v2[jj_r2+s*Nl3] = 3*(jj_r2+s*Nl3) + 2
# T4_full = (4*(3*(jj_r2+s*Nl3)+2) + 2) % N = (12*(jj_r2+s*Nl3) + 10) % N
# This is a v1-element (r=1 mod 3): v1-subindex = (T4_full - 1) / 3 = (12*(jj_r2+s*Nl3)+9)/3 = 4*(jj_r2+s*Nl3)+3
# Reduced mod Nl: (4*jj_r2+3) % Nl for each position.
# sc-type within v1: (4*jj_r2+3) % 3 = (jj_r2+0) % 3 ... for j3%3=2: (4*2+3)%3=11%3=2.
# So Z' for v2-col-r2 comes from sc2 of v1.
# v1[sc2] = (A/rho)*v0[sigma1(sc2)] + 0 (r=1 elements only have the A term)
# sigma1 maps sc2 of v1 to sc1 of v0 (as proved in Script 273).

# So Z' = (A/rho) * v0[sigma1(...)] = t * v0[sigma1(4*jj_r2+3)]
# sigma1(s) for v1 sub-index s: (4*s+2) % Nl
# For s = 4*jj_r2+3: sigma1(s) = (4*(4*jj_r2+3)+2) % Nl = (16*jj_r2+14) % Nl
# sc-type: (16*jj_r2+14) % 3 = (jj_r2+2) % 3. For j3%3=2: (2+2)%3=1. ✓ sc1 of v0.

# For the three elements (s=0,1,2 in the column triplet):
# Full-array T4 for v2[jj_r2+s*Nl3]: gives v1 at (4*(jj_r2+s*Nl3)+3) % Nl.
# Then sigma1 gives v0 at (4*(4*(jj_r2+s*Nl3)+3)+2) % Nl = (16*(jj_r2+s*Nl3)+14) % Nl.

Z0_for_v2r2_sigma1 = np.stack([
    (16*jj_r2 + 14) % Nl,
    (16*(jj_r2+Nl3) + 14) % Nl,
    (16*(jj_r2+2*Nl3) + 14) % Nl,
], axis=1)  # v0 sub-indices

Z0_for_v2r2 = np.stack([
    v0[Z0_for_v2r2_sigma1[:, 0]],
    v0[Z0_for_v2r2_sigma1[:, 1]],
    v0[Z0_for_v2r2_sigma1[:, 2]],
], axis=1) * t_val  # scale by t

# Compute within-col variance of Z-inputs
def within_col_cov2(Z_matrix):
    """Z_matrix: (n, 3). Returns per-column CoV² of Z-inputs."""
    m = Z_matrix.mean(axis=1, keepdims=True)
    v_ = ((Z_matrix - m)**2).mean(axis=1)
    return v_ / (m.squeeze()**2 + 1e-300)

cov2_Z_v0r2 = within_col_cov2(Z0_for_v0r2)  # within-col CoV² of v2-Z for v0-col-r2
cov2_Z_v2r2 = within_col_cov2(Z0_for_v2r2)  # within-col CoV² of v0-Z (scaled by t) for v2-col-r2

print("Z-input within-col CoV² for r=2 sub-group:")
print(f"  v0-col-r2 (Z from sc2-v2): mean={cov2_Z_v0r2.mean():.6f}, std={cov2_Z_v0r2.std():.6f}")
print(f"  v2-col-r2 (Z from t*sc1-v0): mean={cov2_Z_v2r2.mean():.6f}, std={cov2_Z_v2r2.std():.6f}")
print(f"  Fraction where Z_v2r2 > Z_v0r2: {np.mean(cov2_Z_v2r2 > cov2_Z_v0r2):.4f}")
print()

# Now check total column CoV² comparison pointwise
col_v0_r2 = np.stack([v0[jj_r2], v0[jj_r2+Nl3], v0[jj_r2+2*Nl3]], axis=1)
col_v2_r2 = np.stack([v2[jj_r2], v2[jj_r2+Nl3], v2[jj_r2+2*Nl3]], axis=1)
cov2_v0_r2 = col_cov2(col_v0_r2)
cov2_v2_r2 = col_cov2(col_v2_r2)
m2m_v0_r2  = col_m2m(col_v0_r2)
m2m_v2_r2  = col_m2m(col_v2_r2)

print("r=2 sub-group POINTWISE comparison:")
print(f"  frac(CoV2_v2>CoV2_v0): {np.mean(cov2_v2_r2 > cov2_v0_r2):.4f}")
print(f"  frac(m2m_v2<m2m_v0):   {np.mean(m2m_v2_r2 < m2m_v0_r2):.4f}")
print(f"  mean(CoV2_v2-CoV2_v0): {np.mean(cov2_v2_r2 - cov2_v0_r2):.6f}")
print(f"  mean(m2m_v0-m2m_v2):   {np.mean(m2m_v0_r2 - m2m_v2_r2):.6f}")
print()

# Correlation of CoV2 diff and m2m diff within r=2
diff_cov2 = cov2_v2_r2 - cov2_v0_r2
diff_m2m  = m2m_v0_r2  - m2m_v2_r2
corr_r2 = float(np.corrcoef(diff_cov2, diff_m2m)[0,1])
print(f"  corr(CoV2_diff, m2m_diff) within r=2: {corr_r2:.4f}")
print()

# ======================================================================
# SECTION 7: SUMMARY TABLE
# ======================================================================
print("="*70)
print("SECTION 7: SUMMARY TABLE — IS THE INEQUALITY STRUCTURAL?")
print()
print(f"{'lam':>5} {'k':>3} {'t':>6} {'t_crit':>7} {'t<t_crit':>9} "
      f"{'frac_CoV2':>10} {'frac_m2m':>9} {'gap':>8}")
for res in results:
    t = res['t']; t_crit = res['t_crit']
    t_ok = t < t_crit
    print(f"{res['lam']:>5.2f} {res['k']:>3d} {t:>6.4f} {t_crit:>7.4f} "
          f"{'YES' if t_ok else 'NO':>9} "
          f"{res['frac_cov2']:>10.4f} {res['frac_m2m']:>9.4f} {res['global_gap']:>8.5f}")

print()
print("CONCLUSION:")
print("  If t < t_crit = sqrt(lam^2-1) for all test cases, the")
print("  GLOBAL Var(v2)>Var(v0) condition holds analytically.")
print("  The POINTWISE fraction shows per-column Jensen applicability.")
print()
print("done")
