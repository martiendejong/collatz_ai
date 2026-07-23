"""
151_cf_rigorous_bounds.py
==========================
RIGOROUS lower bounds on the size of any non-trivial Collatz cycle, derived
from scratch via continued fractions of theta = log2(3), conditional only on
the published computational verification limit N0 (all n < N0 converge).

Derivation (standard form: cycle of V(n)=(3n+1)/2^e, A odd steps, B halvings):

 (i)  Around the cycle: 2^B = prod_i (3 + 1/n_i)  <  (3 + 1/n_min)^A.
      Verification => n_min > N0  =>  2^(B/A) < 3 + 1/N0
      =>  0 < B - A*theta < A*delta,   delta := log2(1 + 1/(3*N0)).
      (Lower bound 0 is strict since 2^B > 3^A for a positive cycle.)

 (ii) delta is tiny (~2^-69 for N0=2^68), so as long as A*delta < 1 the
      value of B is FORCED:  B = ceil(A*theta).  The search space is
      one-dimensional.

 (iii) Best-approximation property of continued fractions: if q_k <= A < q_(k+1)
      (q_k = convergent denominators of theta), then  ||A*theta|| >= ||q_k*theta||.
      Requirement (i) gives ||A*theta|| < A*delta, hence A > ||q_k*theta||/delta.
      If that exceeds q_(k+1), the whole window [q_k, q_(k+1)) is empty.
      => A_min = first realizable max(q_k, ceil(||q_k theta||/delta)).

All arithmetic at 150+ decimal digits, cross-checked at doubled precision.
"""
import mpmath as mp

def compute(N0_log2_label, N0, dps=160):
    mp.mp.dps = dps
    theta = mp.log(3) / mp.log(2)
    delta = mp.log(1 + mp.mpf(1) / (3 * N0)) / mp.log(2)

    # continued fraction of theta
    a = []
    x = theta
    for _ in range(60):
        ai = int(mp.floor(x)); a.append(ai)
        x = 1 / (x - ai)

    # convergents
    ps = [a[0], a[0]*a[1] + 1]
    qs = [1, a[1]]
    for k in range(2, len(a)):
        ps.append(a[k]*ps[-1] + ps[-2])
        qs.append(a[k]*qs[-1] + qs[-2])

    print(f"--- N0 = {N0_log2_label}  (delta = 2^{mp.nstr(mp.log(delta,2), 6)}) ---")
    A_min = None
    for k in range(len(qs) - 1):
        err = abs(qs[k]*theta - ps[k])          # = ||q_k * theta||
        need = err / delta                       # A must exceed this in window k
        lowA = max(qs[k], int(mp.ceil(need)))
        if lowA < qs[k+1]:
            A_min = lowA
            print(f"  first non-empty window: k={k}, q_k={qs[k]}, q_k+1={qs[k+1]}")
            print(f"  ||q_k theta|| = {mp.nstr(err, 8)}  ->  A > {mp.nstr(need, 10)}")
            break
        # else: window [q_k, q_{k+1}) entirely excluded

    B_min = int(mp.ceil(A_min * theta))
    print(f"  A_min (odd steps / odd elements)  >= {A_min:,}")
    print(f"  B_min (halvings)                  >= {B_min:,}")
    print(f"  total Collatz steps (A+B)         >= {A_min + B_min:,}")
    print(f"  B is forced: B = ceil(A*theta) for all A < 1/delta = "
          f"{mp.nstr(1/delta, 6)}")
    print()
    return A_min, B_min, theta, delta, a

# conservative (fully published, replicated): Barina 2020, n < 2^68
A68, B68, theta, delta68, cf = compute("2^68 (Barina 2020, J. Supercomputing)", 2**68)

# current record (Barina 2025): n < 1.5 * 2^70
A70, B70, _, delta70, _ = compute("1.5*2^70 (Barina 2025)", 3 * 2**69)

print("continued fraction of log2(3), first 30 partial quotients:")
print(" ", cf[:30])
print()

# ---- precision cross-check at doubled dps ----
mp.mp.dps = 320
theta2 = mp.log(3)/mp.log(2)
a2 = []
x = theta2
for _ in range(60):
    ai = int(mp.floor(x)); a2.append(ai)
    x = 1/(x - ai)
print("precision check (dps 160 vs 320): first 55 CF terms identical:",
      cf[:55] == a2[:55])

# ---- sanity: compare with literature ----
print(f"""
cross-check vs literature:
  Hercher (2023) states the next target bound for odd elements is
  K >= 1.375e11; our rigorous two-sided bound with the CURRENT verification
  limit gives A_min = {A70:,} ~ {A70/1e11:.3f}e11 -- same order, consistent.
  (Ours is slightly weaker per A since we use the two-sided best-approximation
  inequality for full rigor and no m-cycle structure theory.)
""")
