"""
195_chain_chebyshev.py
======================
TASK 4 of density_one.tex (the along-chain Chebyshev recursion, Lemma D),
made explicit and measured.

Structure (summation-lemma.md, Assembly A'):
  - G(m) = F(4m+2) - F(m), the backbone log-ratio; feed domination at
    threshold eps is EXACTLY the event {G <= -t0(eps)} with
    t0(eps) = -log2(eps * rho * lam^2) > 0, valid for eps < lam^-2/rho.
  - Chebyshev in flow measure W = v / sum v:
       W{G <= -t0} <= Var_W(G) / (t0 + E_W[G])^2 =: delta0(eps)
    (E_W[G] enters with its sign; measured positive = it helps).
  - Chain step: consecutive chain levels are separated by one selected
    feed edge; conditioning on domination at the current level perturbs
    the next level's G-statistics by Law-A-decay factors. The script
    measures those perturbations DIRECTLY:
       mean shift  m1(g) = E_W[G(next) | chain >= g] - E_W[G]
       var ratio   r2(g) = Var_W(G(next) | chain >= g) / Var_W(G)
    and the realized per-step chain ratios r(g) = W{>=g+1}/W{>=g},
    to be compared with delta(eps) = delta0 * (1 + correction).

Output per k: E_W[G], Var_W(G), then per eps: t0, delta0, measured
one-step mass W{G<=-t0}, chain ratios r(1..5), worst conditional
perturbations. PASS criterion (Lemma D): max_g r(g) <= delta0 with the
measured perturbations small (the recursion closes with explicit
constants).
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


for k in (11, 12, 13, 14, 15):
    lam, v = edge_vector(k)
    N, Nl, T4, m1, R1, m3, R3 = make_maps(k)
    A, B1, B3 = lam ** -2.0, lam ** (ALPHA - 2.0), lam ** (ALPHA - 1.0)
    cb = np.minimum(np.minimum(v[:Nl], v[Nl:2 * Nl]), v[2 * Nl:])

    W = v / v.sum()
    F = np.log2(v)
    G = F[T4] - F                       # backbone log-ratio
    EG = float((W * G).sum())
    VG = float((W * (G - EG) ** 2).sum())

    # selected-feed successor (argmin lift of the feed base), as in 188
    stack = np.stack([v[:Nl], v[Nl:2 * Nl], v[2 * Nl:]])
    amin = stack.argmin(axis=0)
    r1full = np.zeros(N, dtype=np.int64)
    r1full[m1] = R1[m1]
    r1full[m3] = R3[m3]
    has = m1 | m3
    tgt = np.full(N, -1, dtype=np.int64)
    tgt[has] = r1full[has] + amin[r1full[has]] * Nl

    print(f"\nk={k}  lam={lam:.6f}  E_W[G]={EG:+.4f}  Var_W(G)={VG:.4f}  "
          f"eps domain < lam^-2 = {lam**-2:.3f}")
    print("    eps    t0     delta0  W{G<=-t0}  r(1)   r(2)   r(3)   "
          "r(4)   r(5)   max|m1(g)|  max r2(g)")
    for eps in (0.05, 0.10, 0.15, 0.20, 0.25):
        t0 = -log2(eps * lam ** 2)
        if t0 + EG <= 0:
            print(f"    {eps:.2f}  t0+E_W[G] <= 0 -- outside Chebyshev domain")
            continue
        delta0 = VG / (t0 + EG) ** 2
        dom = G <= -t0
        mass1 = float(W[dom].sum())
        # chain lengths along the selected-feed path
        glen = dom.astype(np.int64).copy()
        alive = dom & (tgt >= 0)
        pos = np.where(alive, tgt, 0)
        # record, per chain depth g, the set of NEXT-level nodes for the
        # conditional-perturbation measurement
        cond_stats = []
        for depth in range(1, 8):
            # conditional G-statistics at the current chain frontier,
            # weighted by the ORIGIN flow W (the chain's carrying measure)
            if alive.any():
                wsel = W[alive]
                gnext = G[pos[alive]]
                m1g = float((wsel * gnext).sum() / wsel.sum()) - EG
                v1g = float((wsel * (gnext - (m1g + EG)) ** 2).sum()
                            / wsel.sum()) / VG
                cond_stats.append((m1g, v1g))
            nxt_dom = alive & dom[pos]
            glen[nxt_dom] += 1
            nxt_ok = nxt_dom & (tgt[pos] >= 0)
            pos = np.where(nxt_ok, tgt[pos], pos)
            alive = nxt_ok
            if not alive.any():
                break
        flows = [float(W[glen >= j].sum()) for j in range(1, 7)]
        ratios = [flows[j + 1] / flows[j] if flows[j] > 0 else float("nan")
                  for j in range(5)]
        mm = max(abs(m) for m, _ in cond_stats) if cond_stats else 0.0
        rr = max(r for _, r in cond_stats) if cond_stats else 0.0
        ok = all(not (r == r) or r <= delta0 for r in ratios)  # nan-safe
        print(f"    {eps:.2f}  {t0:5.2f}  {delta0:6.3f}  {mass1:9.5f}  "
              + "  ".join(f"{r:5.3f}" if r == r else "  -  " for r in ratios)
              + f"   {mm:8.3f}   {rr:7.3f}"
              + ("" if ok else "   *** ratio > delta0 ***"))
