"""
209_regularity_transfer.py
==========================
Obs 414 -- sub-question (B): regularity transfer.

Does operator-level sigma_W/rho = 0.755 (dead-flat across k=12..17,
Obs 405) imply Perron-vector CV bounded away from zero?

Route: structural row-sum heterogeneity of the K-L operator at fixed lambda.

The K-L operator update (per node i with top-digit residue r_i):
  rho * v_i = A * v_{T4(i)}              [always]
             + B3 * cb_{R3(i)}           [if r_i = 2]
             + B1 * cb_{R1(i)}           [if r_i = 0]

Since B3 > 0 and B1 > 0 and the Perron eigenvector is strictly positive
(Perron-Frobenius on a primitive matrix), the "bonus" terms are
strictly positive for r=0 and r=2 nodes but ZERO for r=1 nodes.

Structural argument:
  type-2 bonus = B3 * cb_mean > 0
  type-1 bonus = 0
  Delta_r = type-2 bonus - type-1 bonus = B3 * cb_mean > 0

This heterogeneity is a function of lambda ALONE -- independent of k.
Therefore:
  CV(v^(k)) >= c * Delta_r / rho  for some k-independent c > 0.

Measurement plan:
  (A) Bonus spread: at each k, compute mean bonus per node type,
      compare to rho. This is lambda-dependent, k-invariant.
  (B) Cross-lambda grid: sigma_W/rho and CV both measured at
      lambda = 1.30, 1.50, 1.70, 1.90 for k = 14.
      Ratio CV / sigma_W should be stable -> sigma_W is a valid
      proxy for CV's floor.
  (C) Spectral gap: measure (rho - |second eigenvalue|) / rho.
      Bounded gap -> eigenvector nonuniformity is stable.
"""
import numpy as np
from math import log2

ALPHA = log2(3.0)


# ------------------------------------------------------------------ build / perron

def make_params(lam):
    return lam ** -2.0, lam ** (ALPHA - 2.0), lam ** (ALPHA - 1.0)


def build(k, lam):
    A, B1, B3 = make_params(lam)
    N = 3 ** (k - 1)
    i = np.arange(N, dtype=np.int64)
    T4 = (4 * i + 2) % N
    s, r = np.divmod(i, 3)
    Nl = N // 3
    m2 = (r == 2)
    m0 = (r == 0)
    R3 = (2 * s + 1) % Nl   # type-2 lift target
    R1 = (4 * s) % Nl       # type-0 lift target
    return N, Nl, T4, m0, m2, R1, R3, A, B1, B3


def perron(k, lam, n_iter=300):
    N, Nl, T4, m0, m2, R1, R3, A, B1, B3 = build(k, lam)
    v = np.ones(N, dtype=np.float64)
    for _ in range(n_iter):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2 * Nl]), v[2 * Nl:])
        w = A * v[T4]
        w[m2] += B3 * cb[R3[m2]]
        w[m0] += B1 * cb[R1[m0]]
        v = w / w.max()
    # Rayleigh estimate for rho
    cb = np.minimum(np.minimum(v[:Nl], v[Nl:2 * Nl]), v[2 * Nl:])
    w = A * v[T4]
    w[m2] += B3 * cb[R3[m2]]
    w[m0] += B1 * cb[R1[m0]]
    rho = float(w.sum() / v.sum())
    v = v / v.mean()
    return v, rho, cb / cb.mean()


# ------------------------------------------------------------------ Part A: bonus spread

def bonus_spread(k, lam):
    """
    For each node type (0, 1, 2), compute mean 'bonus' in Perron equation.
    Bonus(r=2) = B3 * cb[R3] / rho
    Bonus(r=0) = B1 * cb[R1] / rho
    Bonus(r=1) = 0
    The gap bonus_2 - bonus_1 is the structural source of eigenvector spread.
    """
    N, Nl, T4, m0, m2, R1, R3, A, B1, B3 = build(k, lam)
    v, rho, cb = perron(k, lam)
    # cb is already mean-normalised Nl-vector
    cb_raw = cb * v[:Nl].mean() / cb.mean()  # re-scale: not needed, use raw cb
    # recompute cb from v
    cb2 = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
    bonus2 = float(np.mean(B3 * cb2[R3[m2]])) / rho    # mean bonus for type-2
    bonus0 = float(np.mean(B1 * cb2[R1[m0]])) / rho    # mean bonus for type-0
    # CV of v
    cv = float(np.std(v) / np.mean(v))
    # type-conditional means of v
    Nl = N // 3
    s_all = np.arange(N) % 3  # residue within full vector
    v2_mean = float(np.mean(v[s_all == 2]))
    v1_mean = float(np.mean(v[s_all == 1]))
    v0_mean = float(np.mean(v[s_all == 0]))
    delta_r = bonus2 - 0.0     # gap between type-2 and type-1
    return cv, bonus2, bonus0, delta_r, v2_mean / v1_mean


# ------------------------------------------------------------------ Part B: cross-lam

def sigma_W(k, lam, sel, n_iter=300, seed=1):
    N, Nl, T4, m0, m2, R1, R3, A, B1, B3 = build(k, lam)
    tgt2 = R3 + sel[R3] * Nl
    tgt0 = R1 + sel[R1] * Nl

    def PW(x):
        m = (x[:Nl] + x[Nl:2*Nl] + x[2*Nl:]) / 3.0
        y = x.copy()
        y[:Nl] -= m; y[Nl:2*Nl] -= m; y[2*Nl:] -= m
        return y

    rng = np.random.default_rng(seed)
    d = PW(rng.standard_normal(N))
    d /= np.linalg.norm(d)
    rates = []
    for _ in range(n_iter):
        y = A * d[T4]
        y[m2] += B3 * d[tgt2[m2]]
        y[m0] += B1 * d[tgt0[m0]]
        y = PW(y)
        nrm = np.linalg.norm(y)
        rates.append(nrm)
        d = y / nrm
    tail = np.array(rates[-80:])
    return float(np.exp(np.mean(np.log(tail))))


def cross_lam(lam, k=14):
    v, rho, _ = perron(k, lam)
    cv = float(np.std(v) / np.mean(v))
    # sel for sigma_W
    N = v.size
    Nl = N // 3
    stack = np.stack([v[:Nl], v[Nl:2*Nl], v[2*Nl:]])
    sel = stack.argmin(axis=0).astype(np.int64)
    sw = sigma_W(k, lam, sel)
    _, B1, B3 = make_params(lam)
    delta = (B3 - B1) / rho    # raw spread proxy
    return rho, cv, sw / rho, cv / (sw / rho), delta


# ------------------------------------------------------------------ Part C: spectral gap

def spectral_gap(k, lam, n_iter=200):
    """
    Deflate out the Perron mode; measure the second-largest eigenvalue
    by running the power method on v - <v, ones> * ones (orthogonal
    complement of the Perron direction after normalisation).
    """
    v, rho, _ = perron(k, lam)
    N, Nl, T4, m0, m2, R1, R3, A, B1, B3 = build(k, lam)
    cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])

    def Mv(d):
        w = A * d[T4]
        w[m2] += B3 * cb[R3[m2]] / v[m2] * d[m2]   # linearized
        w[m0] += B1 * cb[R1[m0]] / v[m0] * d[m0]
        return w

    # random start orthogonal to Perron
    rng = np.random.default_rng(42)
    d = rng.standard_normal(N)
    d -= d.dot(v) / v.dot(v) * v
    d /= np.linalg.norm(d)
    rates = []
    for _ in range(n_iter):
        y = Mv(d)
        y -= y.dot(v) / v.dot(v) * v
        nrm = np.linalg.norm(y)
        rates.append(nrm)
        d = y / (nrm + 1e-300)
    tail = np.array(rates[-60:])
    rho2 = float(np.exp(np.mean(np.log(tail + 1e-300))))
    return rho, rho2, (rho - rho2) / rho


# ================================================================== MAIN

LAM = 1.70

print("=" * 68, flush=True)
print("Part A: structural bonus spread vs Perron CV  (lam=1.70)", flush=True)
print("  k   CV(v)    bonus2/rho  bonus0/rho  Delta2/rho  v2/v1", flush=True)
for k in (12, 13, 14, 15, 16):
    cv, b2, b0, dr, ratio = bonus_spread(k, LAM)
    print(f"  {k:2d}  {cv:.5f}  {b2:.5f}     {b0:.5f}     {dr:.5f}     {ratio:.4f}",
          flush=True)

print(flush=True)
print("Part B: cross-lambda calibration  (k=14)", flush=True)
print("  lam   rho       CV(v)    sw/rho   CV/sw    Delta_r/rho", flush=True)
for lam in (1.30, 1.50, 1.70, 1.90):
    rho, cv, sw_r, ratio, delta = cross_lam(lam, k=14)
    print(f"  {lam:.2f}  {rho:.5f}  {cv:.5f}  {sw_r:.5f}  {ratio:.4f}  {delta:.5f}",
          flush=True)

print(flush=True)
print("Part C: spectral gap  (lam=1.70)", flush=True)
print("  k    rho       rho2      gap/rho", flush=True)
for k in (12, 13, 14, 15):
    rho, rho2, gap = spectral_gap(k, LAM)
    print(f"  {k:2d}   {rho:.5f}  {rho2:.5f}  {gap:.5f}", flush=True)

print("done", flush=True)
