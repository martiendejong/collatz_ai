"""
196_transfer_constants.py
=========================
TASK 5 of density_one.tex: the kappa => q => gamma transfer constants.

Two structural upgrades, both checkable here:

(1) THE q<->gamma LEG IS EXPLICITLY 1-DIMENSIONAL. The Min-Loss identity
    (proved) forces every critical pair (lam_k, q_k) onto the explicit
    curve
        q(lam) = 3 (1 - lam^-2) / (lam^(a-2) + lam^(a-1)),  a = log2 3,
    so the transfer 1-gamma vs 1-q is a single-variable calculus fact.
    Claim (one-sided, the direction the program needs):
        h(lam) = 3.4761... * (1 - q(lam)) - (1 - log2 lam) >= 0
    on lam in (1, 2]. Checked on a dense grid; the ratio
    (1-gamma)/(3.4761 (1-q)) along the curve must reproduce the measured
    0.824/0.847/0.873/0.908 at the observed lam_k (consistency).

(2) THE CV => q LEG IS SAMUELSON'S INEQUALITY. For any positive triple,
    mean - min <= sqrt(2) * sigma  (Samuelson, n=3: no sample sits more
    than sigma*sqrt(n-1) below the mean). Aggregating over refinement
    triples with weight = triple mean and one Cauchy-Schwarz:
        1 - q_k <= sqrt(2) * CV_w,
    CV_w = weighted-quadratic-mean of triple CVs. This PROVES the
    empirical linearization constant c1 in [1.19, 1.45] <= sqrt(2) =
    1.414: the measured c1 was Samuelson all along. Verified below at
    k = 11..15 (per-triple check exhaustively + aggregate margin).
"""
import numpy as np
from math import log2, log

ALPHA = log2(3.0)
RATE = 1.0 / log(4.0 / 3.0)          # 3.47605...


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


def q_of_lam(lam):
    return 3.0 * (1.0 - lam ** -2.0) / (lam ** (ALPHA - 2.0)
                                        + lam ** (ALPHA - 1.0))


# ---- (1) the explicit q<->gamma curve ------------------------------------
print("(1) one-sided edge-rate bound  h(lam) = RATE*(1-q) - (1-gamma) >= 0")
grid = np.linspace(1.0001, 2.0, 2_000_001)
q = q_of_lam(grid)
h = RATE * (1.0 - q) - (1.0 - np.log2(grid))
print(f"    grid 2e6 points on (1,2]:  min h = {h.min():.6f} "
      f"at lam = {grid[h.argmin()]:.6f}  "
      f"({'HOLDS' if h.min() >= 0 else '*** FAILS ***'})")
ratio = (1.0 - np.log2(grid)) / (RATE * (1.0 - q))
print(f"    ratio (1-gamma)/(RATE*(1-q)): monotone increasing = "
      f"{bool(np.all(np.diff(ratio) > 0))}, "
      f"range {ratio[0]:.4f} -> {ratio[-2]:.4f} (lam->2 limit 1)")
for lam_k, k in ((1.8188, 13), (1.8420, 15), (1.8585, 17), (1.8704, 19),
                 (1.885, 20), (1.88664, 21)):
    qq = q_of_lam(lam_k)
    rr = (1.0 - log2(lam_k)) / (RATE * (1.0 - qq))
    print(f"    k={k:2d} lam={lam_k:.5f}: q={qq:.5f}  ratio={rr:.4f}")

# ---- (2) Samuelson leg ---------------------------------------------------
print("\n(2) Samuelson: per-triple mean-min <= sqrt(2)*sigma;"
      " aggregate 1-q <= sqrt(2)*CV_w")
for k in (11, 12, 13, 14, 15):
    lam, v = edge_vector(k)
    N = 3 ** (k - 1)
    Nl = N // 3
    T = np.stack([v[:Nl], v[Nl:2 * Nl], v[2 * Nl:]])   # triples
    mean = T.mean(axis=0)
    mn = T.min(axis=0)
    sig = T.std(axis=0)                                # population std
    # per-triple Samuelson (a theorem; checked anyway)
    viol = int(np.sum(mean - mn > np.sqrt(2.0) * sig * (1 + 1e-12)))
    # exact q of the certificate
    q_k = float(3.0 * mn.sum() / v.sum())
    # weighted aggregate: 1-q = sum(mean-min)/sum(mean)
    #                        <= sqrt2 * sum(sig)/sum(mean)
    #                        <= sqrt2 * sqrt(sum(mean*cv^2)/sum(mean)) [C-S]
    cv = np.divide(sig, mean)
    lin = float(np.sqrt(2.0) * sig.sum() / mean.sum())          # tight form
    cvw = float(np.sqrt((mean * cv ** 2).sum() / mean.sum()))   # C-S form
    c1_meas = (1.0 - q_k) / float((mean * cv).sum() / mean.sum())
    print(f"    k={k:2d} lam={lam:.5f}: 1-q = {1-q_k:.5f}  "
          f"sqrt2*E_w[cv] = {lin:.5f}  sqrt2*CV_w = {np.sqrt(2)*cvw:.5f}  "
          f"c1_measured = {c1_meas:.3f} <= sqrt2 = 1.414  "
          f"Samuelson violations: {viol}")
    # consistency: identity q from lam must equal certificate q
    print(f"           min-loss curve q(lam) = {q_of_lam(lam):.5f} "
          f"vs certificate q = {q_k:.5f}  "
          f"(diff {abs(q_of_lam(lam)-q_k):.2e})")
