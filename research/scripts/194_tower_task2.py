"""
194_tower_task2.py
==================
TASK 2 of density_one.tex (Law B profile step), made reproducible.

The Obs 383-388 tower measurements were script-inline (archived in
research/sessions/2026-07-22_26_claude-code-session.jsonl); this script
rebuilds them as a standalone, and adds the Task-2 bound check.

Content:
  (A) Independent re-verification of the BLOCK-EQUATION LEMMA (Obs 383):
        V_p(c) = A * V_p(4c+2 mod 3^p) + b(c) * CB_{p-1}(R(c))
      exactly, for all 2 <= p <= k-1, where V_p = mod-3^p block means of
      the Perron field v, CB_{p-1} = mod-3^(p-1) block means of the min
      field cb, b(c) in {B1, 0, B3} by c mod 3, and
      R(c) = 4*(c div 3) mod 3^(p-1)   (type D1, c = 0 mod 3)
      R(c) = 2*(c div 3)+1 mod 3^(p-1) (type D3, c = 2 mod 3).
  (B) Tower increments Xt_p = log2 V_p - log2 V_{p-1} (exact telescoping
      to F = log2 v): counting-measure variance profile, per-level ratio.
  (C) TASK 2 BOUND: Var(Xt_p) <= C_inj * env^(p-1) with
        C_inj = Var_types( log2 1/(1-phibar_t) )   (k-stable branch means)
        env   = (B1+B3)/3                          (rho = 1 at the edge)
      margins reported per level and per k; any violation printed loudly.
  (D) Task-3 coverage re-check on the same fields: max Xt_p <= log2 3,
      min Xt_p >= -log2(rho*lam^2)  (period-3 block-desert bound).
"""
import numpy as np
from math import log2

ALPHA = log2(3.0)


def make_maps(k):
    N = 3 ** (k - 1)
    i = np.arange(N, dtype=np.int64)
    T4 = (4 * i + 2) % N
    s, r = np.divmod(i, 3)
    Nl = N // 3
    R1 = (4 * s) % Nl
    R3 = (2 * s + 1) % Nl
    return N, Nl, T4, (r == 0), R1, (r == 2), R3


def edge_vector(k):
    N, Nl, T4, m1, R1, m3, R3 = make_maps(k)
    lo_l, hi_l = 1.5, 1.999
    v = np.ones(N, dtype=np.float64)
    for _ in range(36):
        lam = 0.5 * (lo_l + hi_l)
        A, B1, B3 = lam ** -2.0, lam ** (ALPHA - 2.0), lam ** (ALPHA - 1.0)
        w = v.copy()
        for _ in range(60):
            cb = np.minimum(np.minimum(w[:Nl], w[Nl:2 * Nl]), w[2 * Nl:])
            w2 = A * w[T4]
            w2[m1] += B1 * cb[R1[m1]]
            w2[m3] += B3 * cb[R3[m3]]
            g = w2.max()
            w = w2 / g
        if g >= 1.0:
            lo_l, v = lam, w
        else:
            hi_l = lam
    lam = lo_l
    A, B1, B3 = lam ** -2.0, lam ** (ALPHA - 2.0), lam ** (ALPHA - 1.0)
    for _ in range(300):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2 * Nl]), v[2 * Nl:])
        w2 = A * v[T4]
        w2[m1] += B1 * cb[R1[m1]]
        w2[m3] += B3 * cb[R3[m3]]
        v = w2 / w2.max()
    return lam, v


def block_means(field, mod):
    """mod-`mod` block means of `field` (blocks = residue classes)."""
    n = field.size
    idx = np.arange(n, dtype=np.int64) % mod
    sums = np.bincount(idx, weights=field, minlength=mod)
    return sums / (n // mod)


for k in (11, 12, 13, 14, 15):
    lam, v = edge_vector(k)
    N, Nl, T4, m1, R1, m3, R3 = make_maps(k)
    A, B1, B3 = lam ** -2.0, lam ** (ALPHA - 2.0), lam ** (ALPHA - 1.0)
    cb = np.minimum(np.minimum(v[:Nl], v[Nl:2 * Nl]), v[2 * Nl:])

    # ---- (A) block-equation lemma, independent recheck --------------------
    if k == 13:
        print(f"\n(A) block-equation lemma recheck, k={k} (lam={lam:.6f})")
        for p in range(2, k):
            M = 3 ** p
            Vp = block_means(v, M)
            CB = block_means(cb, M // 3)
            c = np.arange(M, dtype=np.int64)
            t = c % 3
            c3 = c // 3
            b = np.where(t == 0, B1, np.where(t == 2, B3, 0.0))
            R = np.where(t == 0, (4 * c3) % (M // 3), (2 * c3 + 1) % (M // 3))
            rhs = A * Vp[(4 * c + 2) % M] + b * CB[R]
            err = np.max(np.abs(rhs - Vp) / Vp)
            print(f"    p={p:2d}: max rel err = {err:.2e}")

    # ---- injection constant from k-stable branch means --------------------
    rho_v = A * v[T4]
    feed = np.zeros(N)
    feed[m1] = B1 * cb[R1[m1]]
    feed[m3] = B3 * cb[R3[m3]]
    phi = feed / (rho_v + feed)
    i3 = np.arange(N, dtype=np.int64) % 3
    phibar_t = [phi[i3 == t].mean() for t in range(3)]
    lifts = [log2(1.0 / (1.0 - pt)) if pt < 1 else float("inf")
             for pt in phibar_t]
    C_inj = float(np.var(lifts))          # population var, uniform 1/3
    env = (B1 + B3) / 3.0
    phibar = float(phi.mean())

    # ---- (B)+(C) tower profile and Task-2 bound ---------------------------
    print(f"\nk={k}  lam={lam:.6f}  env=(B1+B3)/3={env:.4f}  "
          f"phibar={phibar:.4f}  C_inj={C_inj:.4f} bits^2  "
          f"(type lifts {lifts[0]:.3f}/{lifts[1]:.3f}/{lifts[2]:.3f})")
    print("    p   Var(Xt_p)  bound C*env^(p-1)  margin  ratio   "
          "maxXt   minXt")
    Vprev = np.full(1, v.mean())
    var_prev = None
    var_sum = 0.0
    ok = True
    for p in range(1, k):
        M = 3 ** p
        Vp = block_means(v, M)
        Xp = np.log2(Vp) - np.log2(Vprev[np.arange(M) % (M // 3)])
        var_p = float(np.var(Xp))
        var_sum += var_p
        bound = C_inj * env ** (p - 1)
        margin = bound / var_p if var_p > 0 else float("inf")
        ratio = var_p / var_prev if var_prev else float("nan")
        flag = "" if var_p <= bound else "  *** VIOLATION ***"
        if var_p > bound:
            ok = False
        print(f"    {p:2d}  {var_p:9.5f}  {bound:12.5f}     "
              f"{margin:6.2f}  {ratio:5.3f}  {float(Xp.max()):6.3f}  "
              f"{float(Xp.min()):6.3f}{flag}")
        Vprev = Vp
        var_prev = var_p
    F = np.log2(v)
    varF = float(np.var(F))
    neg_bound = log2(lam ** 2)            # rho = 1: -log2(rho*lam^2)
    print(f"    Var(F)={varF:.4f}  sum Var(Xt_p)={var_sum:.4f}  "
          f"(cross-terms {'NEGATIVE (sum >= VarF)' if var_sum >= varF else 'POSITIVE'})")
    print(f"    Task-3 coverage: max Xt <= log2 3 = {ALPHA:.3f}; "
          f"min Xt >= -log2(lam^2) = -{neg_bound:.3f}")
    print(f"    TASK 2 BOUND: {'HOLDS at every level' if ok else 'FAILS'}")
