"""
162_krasikov_system.py
=======================
CAMPAIGN STEP 1: an operational, rigorously-conservative Krasikov-type
difference-inequality system for lower bounds f(x) >> x^gamma on the count
of integers reaching 1.

Scheme (conservative variant of Applegate-Lagarias / Krasikov-Lagarias):
  - States: residue classes mod 3^tau.
  - For source class s, enumerate ALL 3^j lifts to mod 3^(tau+j); from each
    lift (as an exact integer representative) enumerate all inverse-tree
    paths of length j (children n = (m*2^e-1)/3, e <= E, parity forced by
    m mod 3, interior nodes != 0 mod 3), recording leaf class mod 3^tau and
    size factor prod(2^e_i/3).
  - RIGOR SAFEGUARDS (each one only weakens the bound, never breaks it):
      (a) entrywise MIN over lifts -> a single matrix valid for every
          integer in class s;
      (b) prune paths with factor <= 1 (keeps the downward induction sound);
      (c) e-cap E discards mass.
  - Certificate: gamma* = sup{gamma : spectral radius rho(M(gamma)) >= 1},
    where M(gamma)[c,s] = min over lifts of sum over paths factor^(-gamma).
    Then f_c(x) >= const * x^gamma* by standard induction.

Historic ladder for calibration: Krasikov 1989: 0.43; tree-search 1995:
0.654; Krasikov inequalities mod 3^9 (LP, 1995): 0.809; nonlinear 2003: 0.84.
Our min-over-lifts crush is strictly weaker than their LP -- the goal here is
an operational, scaling baseline, not yet parity with 0.84.
"""
import numpy as np
from itertools import product
from math import log2

def build_matrix_paths(tau, j, E):
    """For each source class s mod 3^tau: min over lifts of the path-mass
    vector into target classes. Returns list of dicts or the raw per-lift
    data; we store per (s, lift) a list of (target_class, log2_factor)."""
    T3 = 3**tau
    L3 = 3**j
    LOG23 = log2(3)
    per_source = []
    for s in range(T3):
        lift_results = []
        for t in range(L3):
            m0 = s + T3 * t
            if m0 % 3 == 0:
                # a root == 0 mod 3 has no children; contributes nothing
                lift_results.append([])
                continue
            if m0 == 0:
                lift_results.append([])
                continue
            # DFS: (value, depth, log2factor)
            out = []
            stack = [(m0, 0, 0.0)]
            while stack:
                m, d, lf = stack.pop()
                if d == j:
                    out.append((m % T3, lf))
                    continue
                r = m % 3
                if r == 0:
                    continue  # interior must branch; dead path
                e0 = 2 if r == 1 else 1
                for e in range(e0, E + 1, 2):
                    n = (m * (1 << e) - 1) // 3
                    stack.append((n, d + 1, lf + e - LOG23))
            lift_results.append(out)
        per_source.append(lift_results)
    return per_source

def matrix_at_gamma(per_source, tau, gamma):
    T3 = 3**tau
    M = np.zeros((T3, T3))
    for s, lift_results in enumerate(per_source):
        acc = None
        for out in lift_results:
            v = np.zeros(T3)
            for c, lf in out:
                if lf <= 0:          # factor <= 1: prune (induction safety)
                    continue
                v[c] += 2.0 ** (-gamma * lf)
            acc = v if acc is None else np.minimum(acc, v)
        if acc is not None:
            M[:, s] = acc
    return M

def spectral_radius(M, iters=300):
    n = M.shape[0]
    v = np.ones(n) / n
    lam = 0.0
    for _ in range(iters):
        w = M @ v
        nrm = np.linalg.norm(w)
        if nrm == 0:
            return 0.0
        v = w / nrm
        lam = nrm
    return float(lam)

def certificate(tau, j, E):
    per_source = build_matrix_paths(tau, j, E)
    lo, hi = 0.0, 1.0
    # rho decreasing in gamma; find gamma with rho = 1
    for _ in range(40):
        mid = (lo + hi) / 2
        rho = spectral_radius(matrix_at_gamma(per_source, tau, mid))
        if rho >= 1.0:
            lo = mid
        else:
            hi = mid
    return lo

print("Krasikov-type conservative certificates  f(x) >> x^gamma")
print(f"{'tau':>4} {'j':>3} {'E':>3} {'states':>7} {'gamma*':>8}")
for tau, j, E in [(1, 3, 6), (2, 3, 6), (2, 4, 8), (3, 4, 8),
                  (3, 5, 8), (3, 5, 10), (3, 6, 10)]:
    g = certificate(tau, j, E)
    print(f"{tau:>4} {j:>3} {E:>3} {3**tau:>7} {g:>8.4f}", flush=True)

print("""
calibration: tree-search 1995 = 0.654; Krasikov LP mod 3^9 = 0.809;
nonlinear 2003 = 0.84. The min-over-lifts crush costs precision by design;
the scaling of gamma* with (tau, j, E) maps what the FULL inequality system
(replacing min by the Applegate-Lagarias LP over all lift inequalities)
would need to reproduce 0.809 and then attack 0.84 at mod 3^k, k > 9.""")
