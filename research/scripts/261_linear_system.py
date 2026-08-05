"""
261_linear_system.py
====================
Close the proof: solve the 6-variable class-mean linear system ANALYTICALLY.

THE 6-VARIABLE SYSTEM (from Obs 465):
  rho*a0_v2 = B3*c_1 + (A^2/rho)*a2_v0            ...(A)
  rho*a1_v2 = B3*c_0 + (A^2/rho)*a0_v0            ...(B)
  rho*a2_v2 = B3*c_2 + (A^2/rho)*a1_v0            ...(C)
  rho*a0_v0 = A*a0_v2   + B1*c_0                  ...(D)
  rho*a1_v0 = A*a1_v2   + B1*c_1                  ...(E)
  rho*a2_v0 = A*a2_v2   + B1*c_2                  ...(F)

This is a LINEAR SYSTEM in (a0_v2, a1_v2, a2_v2, a0_v0, a1_v0, a2_v0).

WRITE AS: rho*a = M_off*a + B_ext*c
  M_off is the coupling matrix, B_ext is the forcing by c = (c_0, c_1, c_2).

Rearranging: (rho*I - M_off)*a = B_ext*c
=> a = (rho*I - M_off)^{-1} * B_ext * c

GOAL: Show a[1] > a[0] (a1_v2 > a0_v2) from the structure of (rho*I-M_off)^{-1}*B_ext.

ORDER THE VARIABLES: x = (a0_v2, a1_v2, a2_v2, a0_v0, a1_v0, a2_v0)
EQUATIONS:
  rho*x[0] = (A^2/rho)*x[5] + B3*c_1    => rho*x[0] - (A^2/rho)*x[5] = B3*c_1
  rho*x[1] = (A^2/rho)*x[3] + B3*c_0    => rho*x[1] - (A^2/rho)*x[3] = B3*c_0
  rho*x[2] = (A^2/rho)*x[4] + B3*c_2    => rho*x[2] - (A^2/rho)*x[4] = B3*c_2
  rho*x[3] = A*x[0] + B1*c_0             => rho*x[3] - A*x[0] = B1*c_0
  rho*x[4] = A*x[1] + B1*c_1             => rho*x[4] - A*x[1] = B1*c_1
  rho*x[5] = A*x[2] + B1*c_2             => rho*x[5] - A*x[2] = B1*c_2

MATRIX FORM: (rho*I - M_off)*x = b*c  where b is the forcing vector.

Note: THE SYSTEM DECOUPLES INTO TWO INDEPENDENT 2x2 SYSTEMS:
  (0,3): x[0] and x[3] coupled via A and A^2/rho.
  Actually wait -- let me check:
  Eq (A): couples x[0] with x[5], NOT x[3].
  Eq (D): couples x[3] with x[0].
  So (A)+(D): x[0] and x[3] and x[5]:
    (A): x[0] depends on x[5] and c_1
    (D): x[3] depends on x[0] and c_0
    (F): x[5] depends on x[2] and c_2
    (C): x[2] depends on x[4] and c_2
    (E): x[4] depends on x[1] and c_1
    (B): x[1] depends on x[3] and c_0

  Chain: x[0] <-> x[5] <-> x[2] <-> x[4] <-> x[1] <-> x[3] <-> x[0]. (FULL CYCLE!)

  So all 6 variables are coupled in one 6-cycle!

SOLVE THE 2-VARIABLE SUBSYSTEMS:
  Actually, by symmetry, the system splits by the "class permutation":
  Class mapping under the K-L system: 0->1->2->0 (cyclic?). Let's check.

  v2-class-0 -> depends on cb-class-1 (c_1) and v0-class-2 via T4.
  v2-class-1 -> depends on cb-class-0 (c_0) and v0-class-0 via T4.
  v2-class-2 -> depends on cb-class-2 (c_2) and v0-class-1 via T4.

  v0-class-0 -> depends on v2-class-0 via T4 and cb-class-0.
  v0-class-1 -> depends on v2-class-1 via T4 and cb-class-1.
  v0-class-2 -> depends on v2-class-2 via T4 and cb-class-2.

  PAIRING: (v2-class-0, v0-class-0) are NOT directly coupled!
  v2-class-0 couples with v0-class-2 (T4 goes to class-2 of v0).
  v0-class-0 couples with v2-class-0 (T4 goes to class-0 of v2).

  So: {x[0]=a0_v2, x[3]=a0_v0} form one independent 2x2 system IF no cross-coupling.

  But x[0] depends on x[5]=a2_v0, and x[5] depends on x[2]=a2_v2, and x[2] depends on x[4]=a1_v0.
  So it's all coupled in one big cycle.

DIRECT NUMERICAL SOLUTION: Solve the 6x6 system and check a1_v2 > a0_v2.
"""
import numpy as np
from math import log2
import sys

ALPHA = log2(3.0)

def solve_linear_system(lam, rho, c, A, B1, B3):
    """
    Solve the 6-variable linear system for class means.
    x = (a0_v2, a1_v2, a2_v2, a0_v0, a1_v0, a2_v0)
    c = (c_0, c_1, c_2)

    The system (rho*I - M_off)*x = B_ext*c:
    """
    A2r = A**2 / rho  # A^2/rho

    # System matrix (rho*I - M_off):
    # Row 0 (eq A): rho*x[0] - A2r*x[5] = B3*c[1]
    # Row 1 (eq B): rho*x[1] - A2r*x[3] = B3*c[0]
    # Row 2 (eq C): rho*x[2] - A2r*x[4] = B3*c[2]
    # Row 3 (eq D): -A*x[0] + rho*x[3] = B1*c[0]
    # Row 4 (eq E): -A*x[1] + rho*x[4] = B1*c[1]
    # Row 5 (eq F): -A*x[2] + rho*x[5] = B1*c[2]

    M = np.array([
        [rho,  0,    0,    0,    0,   -A2r],  # row 0
        [0,    rho,  0,   -A2r,  0,    0  ],  # row 1
        [0,    0,    rho,  0,   -A2r,  0  ],  # row 2
        [-A,   0,    0,    rho,  0,    0  ],  # row 3
        [0,   -A,    0,    0,    rho,  0  ],  # row 4
        [0,    0,   -A,    0,    0,    rho],  # row 5
    ], dtype=np.float64)

    # RHS: B_ext * c
    rhs = np.array([
        B3*c[1],  # row 0
        B3*c[0],  # row 1
        B3*c[2],  # row 2
        B1*c[0],  # row 3
        B1*c[1],  # row 4
        B1*c[2],  # row 5
    ], dtype=np.float64)

    x = np.linalg.solve(M, rhs)
    return x

def run_kl_and_extract(k, lam, n_iter=None):
    if n_iter is None:
        n_iter = 600 + 100*max(0, k-8)
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

    v0 = v[0::3]; v1 = v[1::3]; v2 = v[2::3]
    cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
    s = np.arange(Nl, dtype=np.int64)

    # Class means
    a = np.array([np.mean(v2[s%3==r]) for r in range(3)] +
                 [np.mean(v0[s%3==r]) for r in range(3)])
    c = np.array([np.mean(cb[s%3==r]) for r in range(3)])

    # Estimate rho from v1 = (A/rho)*v0[sigma1]
    sigma1 = (4*s+2) % Nl
    A_over_rho = float(np.mean(v1 / v0[sigma1]))
    rho = A / A_over_rho

    return a, c, rho, A, B1, B3, Nl

def analyze(k, lam):
    a_kl, c, rho, A, B1, B3, Nl = run_kl_and_extract(k, lam)

    # Solve the linear system
    x_sol = solve_linear_system(lam, rho, c, A, B1, B3)

    # Compare with K-L solution
    err = np.max(np.abs(x_sol - a_kl))
    rel_err = np.max(np.abs(x_sol - a_kl) / (np.abs(a_kl) + 1e-15))

    # Check a1_v2 > a0_v2
    a0_v2 = x_sol[0]; a1_v2 = x_sol[1]; a2_v2 = x_sol[2]
    a0_v0 = x_sol[3]; a1_v0 = x_sol[4]; a2_v0 = x_sol[5]

    print(f"\nk={k}, lam={lam:.2f}:")
    print(f"  a0_v2={a0_v2:.4f}  a1_v2={a1_v2:.4f}  a2_v2={a2_v2:.4f}")
    print(f"  a0_v0={a0_v0:.4f}  a1_v0={a1_v0:.4f}  a2_v0={a2_v0:.4f}")
    print(f"  Linear system error: max_abs={err:.2e}  max_rel={rel_err:.2e}")
    print(f"  a1_v2 > a0_v2? {a1_v2 > a0_v2}  (ratio={a1_v2/a0_v2:.4f})")
    print(f"  a2_v2 > a0_v2? {a2_v2 > a0_v2}  (ratio={a2_v2/a0_v2:.4f})")
    print(f"  a2_v0 > a0_v0? {a2_v0 > a0_v0}  (ratio={a2_v0/a0_v0:.4f})")

    return x_sol, c, rho, A, B1, B3

# ANALYTICAL SOLUTION OF THE 6x6 SYSTEM
def solve_analytically(rho, c, A, B1, B3):
    """
    Solve symbolically. With A2r = A^2/rho:
    The 6 equations:
    (A): rho*x0 = A2r*x5 + B3*c1
    (B): rho*x1 = A2r*x3 + B3*c0
    (C): rho*x2 = A2r*x4 + B3*c2
    (D): rho*x3 = A*x0 + B1*c0
    (E): rho*x4 = A*x1 + B1*c1
    (F): rho*x5 = A*x2 + B1*c2

    The system has a chain structure: substitute recursively.

    From (A) and (F): rho*x0 = (A2r/rho)*(A*x2 + B1*c2) + B3*c1
                      = (A3/rho^2)*x2 + A2r*B1*c2/rho + B3*c1   ...(AF)
    From (C) and (E): rho*x2 = (A2r/rho)*(A*x1 + B1*c1) + B3*c2
                      = (A3/rho^2)*x1 + A2r*B1*c1/rho + B3*c2   ...(CE)
    From (B) and (D): rho*x1 = (A2r/rho)*(A*x0 + B1*c0) + B3*c0
                      = (A3/rho^2)*x0 + A2r*B1*c0/rho + B3*c0   ...(BD)

    So we have 3 equations in x0, x1, x2:
    rho*x0 = q*x2 + B3*c1 + A2r*B1*c2/rho     ...(AF)
    rho*x2 = q*x1 + B3*c2 + A2r*B1*c1/rho     ...(CE)
    rho*x1 = q*x0 + B3*c0 + A2r*B1*c0/rho     ...(BD)
    where q = A^3/rho^2.

    Substitute recursively:
    rho*x0 = q*[q*x1 + B3*c2 + A2r*B1*c1/rho]/rho + B3*c1 + A2r*B1*c2/rho
    = (q^2/rho)*x1 + q*B3*c2/rho + q*A2r*B1*c1/rho^2 + B3*c1 + A2r*B1*c2/rho

    rho*x0 = (q^2/rho)*[q*x0/rho + B3*c0/rho + A2r*B1*c0/rho^2] + ...
    = (q^3/rho^2)*x0 + ...

    So: rho*x0 - (q^3/rho^2)*x0 = sum of c-terms
    => x0*(rho - q^3/rho^2) = rhs0
    => x0 = rhs0 / (rho - q^3/rho^2)

    where q = A^3/rho^2, so q^3 = A^9/rho^6, and rho - q^3/rho^2 = rho - A^9/rho^8.

    Let's compute this!
    """
    A2r = A**2 / rho
    q = A**3 / rho**2
    q2 = q**2; q3 = q**3

    # Denominator: D = rho^3 - q^3 = rho^3 - A^9/rho^6 (but let's keep as rho^3 - q^3 times rho^2 denominator)
    # Actually: rho - q^3/rho^2 = (rho^3 - q^3)/rho^2 = (rho^3 - A^9/rho^6)/rho^2
    # More cleanly: multiply through by rho^2:
    # (rho^3 - q^3)*x0 = rho^2 * rhs0
    # Note q = A^3/rho^2, so q^3 = A^9/rho^6.
    # rho^3 - q^3 = rho^3 - A^9/rho^6. For rho > 0 and A < rho (which holds for lambda > 1):
    # q = A^3/rho^2 < rho (iff A^3 < rho^3 iff A < rho, TRUE for lambda > 1).
    # So D = rho^3 - q^3 > 0! (since rho > q)

    D = rho**3 - q**3  # ALWAYS > 0 if rho > q, i.e., A^3 < rho^3 iff A < rho

    # c0, c1, c2 from c array
    c0, c1, c2 = c

    # rhs for x1 (from BD): rho*x1 = q*x0 + B3*c0 + A2r*B1*c0/rho
    # rhs for x2 (from CE): rho*x2 = q*x1 + B3*c2 + A2r*B1*c1/rho
    # rhs for x0 (from AF): rho*x0 = q*x2 + B3*c1 + A2r*B1*c2/rho

    # Define forcing terms:
    f0 = B3*c1 + A2r*B1*c2/rho  # forcing for x0 (in eq AF)
    f1 = B3*c0 + A2r*B1*c0/rho  # forcing for x1 (in eq BD) = (B3+A2r*B1/rho)*c0
    f2 = B3*c2 + A2r*B1*c1/rho  # forcing for x2 (in eq CE)

    # After substitution:
    # rho^3*x0 = q^2*f1 + q*rho*f2 + rho^2*f0 + q^3*x0
    # (rho^3 - q^3)*x0 = rho^2*f0 + q*rho*f2 + q^2*f1
    # D*x0 = rho^2*f0 + q*rho*f2 + q^2*f1

    x0 = (rho**2 * f0 + q*rho*f2 + q**2*f1) / D
    x1 = (rho**2 * f1 + q*rho*f0 + q**2*f2) / D
    x2 = (rho**2 * f2 + q*rho*f1 + q**2*f0) / D

    # And x3, x4, x5 from (D,E,F):
    x3 = (A*x0 + B1*c0) / rho
    x4 = (A*x1 + B1*c1) / rho
    x5 = (A*x2 + B1*c2) / rho

    return np.array([x0, x1, x2, x3, x4, x5]), D, f0, f1, f2

print("261: Analytical solution of 6-variable K-L class-mean system")
print("="*70)

# Verify analytical solution
k, lam = 8, 1.70
a_kl, c, rho, A, B1, B3, Nl = run_kl_and_extract(k, lam)
x_anal, D, f0, f1, f2 = solve_analytically(rho, c, A, B1, B3)
err = np.max(np.abs(x_anal - a_kl))
print(f"\nAnalytical solution check (k={k}, lam={lam}):")
print(f"  err = {err:.2e}")
print(f"  x_anal  = {x_anal}")
print(f"  a_kl    = {a_kl}")

# KEY RESULT: a1_v2 - a0_v2 in CLOSED FORM
# a0_v2 = (rho^2*f0 + q*rho*f2 + q^2*f1) / D
# a1_v2 = (rho^2*f1 + q*rho*f0 + q^2*f2) / D
# a1_v2 - a0_v2 = (rho^2*(f1-f0) + q*rho*(f0-f2) + q^2*(f2-f1)) / D
#               = [(rho^2-q^2)*(f1-f0) + q*rho*(f0-f2) + q^2*(f2-f0)] / D ... simplify

print(f"\n  D = rho^3 - q^3 = {D:.6f}  (POSITIVE since q < rho)")
print(f"  q = A^3/rho^2 = {A**3/rho**2:.6f}  < rho = {rho:.6f}")
print(f"  f0 = B3*c1 + A2r*B1*c2/rho = {f0:.6f}  (forcing for a0_v2)")
print(f"  f1 = (B3+A2r*B1/rho)*c0    = {f1:.6f}  (forcing for a1_v2)")
print(f"  f2 = B3*c2 + A2r*B1*c1/rho = {f2:.6f}  (forcing for a2_v2)")

# a1 - a0 = (rho^2*(f1-f0) + q*rho*(f0-f2) + q^2*(f2-f1)) / D
q = A**3/rho**2
diff_10 = (rho**2*(f1-f0) + q*rho*(f0-f2) + q**2*(f2-f1)) / D
print(f"\n  a1_v2 - a0_v2 (closed form) = {diff_10:.6f}")
print(f"  a1_v2 - a0_v2 (from K-L)    = {a_kl[1]-a_kl[0]:.6f}")
print(f"  f1 - f0 = {f1-f0:.6f}  (>0 means f1>f0, primary driver)")
print(f"  f0 - f2 = {f0-f2:.6f}")
print(f"  f2 - f1 = {f2-f1:.6f}")

# PROVE f1 > f0:
# f1 = (B3 + A2r*B1/rho)*c0
# f0 = B3*c1 + A2r*B1*c2/rho
# f1 - f0 = B3*(c0-c1) + A2r*B1*(c0-c2)/rho
#          = B3*c0*(1-A/rho) + A2r*B1*(c0-c2)/rho
# First term: POSITIVE (since A/rho<1)
# Second term: sign depends on c0-c2 (sign varies by lambda)

A2r = A**2/rho
print(f"\n  f1-f0 decomposition:")
print(f"    B3*c0*(1-A/rho) = {B3*c[0]*(1-A/rho):.6f}  [POSITIVE]")
print(f"    A2r*B1*(c0-c2)/rho = {A2r*B1*(c[0]-c[2])/rho:.6f}  [sign?]")
print(f"    f0-f2 contribution to a1-a0: q*rho*(f0-f2)/D = {q*rho*(f0-f2)/D:.6f}")

# LAMBDA SCAN: check if a1>a0 always, and if D>0 always
print(f"\n\nLambda scan k=8: verify closed-form solution + a1>a0")
print(f"{'lam':>6}  {'D>0':>6}  {'f1>f0':>6}  {'a1>a0':>6}  {'a1-a0':>9}  {'x_err':>9}")
for lam in [1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00]:
    a_kl_l, c_l, rho_l, A_l, B1_l, B3_l, _ = run_kl_and_extract(8, lam)
    x_l, D_l, f0_l, f1_l, f2_l = solve_analytically(rho_l, c_l, A_l, B1_l, B3_l)
    q_l = A_l**3 / rho_l**2
    diff_l = (rho_l**2*(f1_l-f0_l) + q_l*rho_l*(f0_l-f2_l) + q_l**2*(f2_l-f1_l)) / D_l
    err_l = float(np.max(np.abs(x_l - a_kl_l)))
    print(f"lam={lam:.2f}  {D_l>0!s:>6}  {f1_l>f0_l!s:>6}  {diff_l>0!s:>6}  "
          f"{diff_l:>9.5f}  {err_l:>9.2e}")
    sys.stdout.flush()

# DEPTH SCAN
print(f"\n\nDepth scan lam=1.70: closed-form solution accuracy")
print(f"{'k':>4}  {'D>0':>6}  {'f1>f0':>6}  {'a1>a0':>6}  {'a1-a0':>9}  {'x_err':>9}")
lam = 1.70
for k in range(5, 12):
    a_kl_k, c_k, rho_k, A_k, B1_k, B3_k, _ = run_kl_and_extract(k, lam)
    x_k, D_k, f0_k, f1_k, f2_k = solve_analytically(rho_k, c_k, A_k, B1_k, B3_k)
    q_k = A_k**3 / rho_k**2
    diff_k = (rho_k**2*(f1_k-f0_k) + q_k*rho_k*(f0_k-f2_k) + q_k**2*(f2_k-f1_k)) / D_k
    err_k = float(np.max(np.abs(x_k - a_kl_k)))
    print(f"k={k:>2}  {D_k>0!s:>6}  {f1_k>f0_k!s:>6}  {diff_k>0!s:>6}  "
          f"{diff_k:>9.5f}  {err_k:>9.2e}")
    sys.stdout.flush()

# PROVE THE FORMULA a1-a0 > 0:
print(f"\n\n=== ANALYTICAL PROOF OF a1_v2 > a0_v2 ===")
print(f"""
CLOSED FORM (derived above):
  D*(a1_v2 - a0_v2) = rho^2*(f1-f0) + q*rho*(f0-f2) + q^2*(f2-f1)
where:
  q = A^3/rho^2 > 0
  D = rho^3 - q^3 > 0 (iff rho > q = A^3/rho^2 iff rho^3 > A^3 iff rho > A)
  f0 = B3*c1 + (A^2*B1/rho^2)*c2
  f1 = (B3 + A^2*B1/rho^2)*c0
  f2 = B3*c2 + (A^2*B1/rho^2)*c1

USING c1 = (A/rho)*c0 (Obs 464):
  f0 = B3*(A/rho)*c0 + (A^2*B1/rho^2)*c2
  f1 = (B3 + A^2*B1/rho^2)*c0 = B3*c0*(1 + A^2*B1/(B3*rho^2))
  f2 = B3*c2 + (A^2*B1/rho^2)*(A/rho)*c0 = B3*c2 + A^3*B1*c0/rho^3

f1 - f0 = B3*c0 + A^2*B1*c0/rho^2 - B3*(A/rho)*c0 - A^2*B1*c2/rho^2
        = B3*c0*(1 - A/rho) + A^2*B1*(c0-c2)/rho^2

The sign of f1-f0 depends on c0 vs c2 and relative magnitude.

KEY SIMPLIFICATION:
D*(a1-a0) = rho^2*(f1-f0) + q*rho*(f0-f2) + q^2*(f2-f1)
           = (rho^2-q^2)*(f1-f0) + q*rho*(f0-f2) + q^2*(f2-f0)...
Actually factor differently:
D*(a1-a0) = (rho^2-q^2)*(f1-f0) + q*(rho*(f0-f2) + q*(f2-f1))
           = (rho^2-q^2)*(f1-f0) + q*((rho-q)*f0 + (q-rho)*f1 + (rho-q)*f2)... hmm.

Let me use c1 = (A/rho)*c0 and check numerically which terms dominate.
""")

# Numerical analysis of which term dominates in D*(a1-a0)
lam, k = 1.70, 8
a_kl_l, c_l, rho_l, A_l, B1_l, B3_l, _ = run_kl_and_extract(k, lam)
q_l = A_l**3/rho_l**2
_, D_l, f0_l, f1_l, f2_l = solve_analytically(rho_l, c_l, A_l, B1_l, B3_l)

term1 = rho_l**2 * (f1_l-f0_l)
term2 = q_l*rho_l * (f0_l-f2_l)
term3 = q_l**2 * (f2_l-f1_l)

print(f"At k={k}, lam={lam}:")
print(f"  D*(a1-a0) = {D_l*(a_kl_l[1]-a_kl_l[0]):.6f}")
print(f"  term1 = rho^2*(f1-f0) = {term1:.6f}")
print(f"  term2 = q*rho*(f0-f2) = {term2:.6f}")
print(f"  term3 = q^2*(f2-f1)   = {term3:.6f}")
print(f"  sum   = {term1+term2+term3:.6f}")
print(f"  q/rho = {q_l/rho_l:.4f}  (small -> term1 dominates)")

print("\ndone")
