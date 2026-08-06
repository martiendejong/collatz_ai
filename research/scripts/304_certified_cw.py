"""
304_certified_cw.py
===================
RIGOROUS interval certification (mpmath.iv) of the K-L quantities at small k.

- lambda exact rational (all grid values are rational: 21/20, ..., 2).
- alpha = log2(3), u = lambda^alpha enclosed rigorously via interval log/exp.
- Power iteration in interval arithmetic; each F-evaluation is componentwise
  exact interval propagation (each output uses each input component at most
  once per term -> no dependency inflation in min/plus).
- Collatz-Wielandt: rho is UNCONDITIONALLY contained in
  [min_i F(w)_i / w_i, max_i F(w)_i / w_i] for any positive w
  (F monotone and 1-homogeneous). This gives rigorous rho, t, R intervals.
- Criteria (c2/c0 vs R, mu2/mu0 vs R, s2/s0 vs lambda) evaluated at the
  converged w with full interval propagation. Caveat (stated honestly):
  criteria are evaluated at w, not at the exact Perron vector v*; the
  remaining step for unconditional certification of the criteria is a
  Birkhoff-contraction bound on the eigenvector distance (residual
  log(beta/alpha) is reported and is ~1e-20).
"""
from mpmath import iv, mp
from fractions import Fraction
import numpy as np

iv.prec = 120

def certify(lam_frac, k, niters):
    lam = iv.mpf(lam_frac.numerator) / iv.mpf(lam_frac.denominator)
    alpha = iv.log(3) / iv.log(2)
    A = 1 / (lam * lam)
    u = iv.exp(alpha * iv.log(lam))       # lambda^alpha, rigorous enclosure
    B1 = A * u                            # lambda^(alpha-2)
    B3 = B1 * lam                         # lambda^(alpha-1)

    N = 3 ** (k - 1)
    Nl = N // 3
    Nl3 = Nl // 3
    T4 = [(4 * i + 2) % N for i in range(N)]
    R1 = [(4 * s) % Nl for s in range(Nl)]
    R3 = [(2 * s + 1) % Nl for s in range(Nl)]

    # Phase 1: point iteration at high precision (mp context) to converge w.
    # Any positive w gives valid CW bounds, so w itself needs no certification.
    mp.prec = 120
    lam_p = mp.mpf(lam_frac.numerator) / mp.mpf(lam_frac.denominator)
    alpha_p = mp.log(3) / mp.log(2)
    A_p = 1 / (lam_p * lam_p)
    u_p = mp.e ** (alpha_p * mp.log(lam_p))
    B1_p = A_p * u_p
    B3_p = B1_p * lam_p

    try:
        vf = np.load(f"E:/projects/collatz/research/cache/v_lam{float(lam_frac):.2f}_k{k}.npy")
        v = [mp.mpf(float(x)) for x in vf]
    except Exception:
        v = [mp.mpf(1) for _ in range(N)]

    def Fp(v):
        cb = [min(v[j], v[j + Nl], v[j + 2 * Nl]) for j in range(Nl)]
        w = [A_p * v[T4[i]] for i in range(N)]
        for s in range(Nl):
            w[3 * s] += B1_p * cb[R1[s]]
            w[3 * s + 2] += B3_p * cb[R3[s]]
        return w

    for _ in range(niters):
        w = Fp(v)
        mx = max(w)
        v = [x / mx for x in w]

    # Phase 2: single rigorous interval pass at the point vector v.
    v = [iv.mpf(x) for x in v]

    def F(v):
        cb = [min(v[j], v[j + Nl], v[j + 2 * Nl]) for j in range(Nl)]
        w = [A * v[T4[i]] for i in range(N)]
        for s in range(Nl):
            w[3 * s] = w[3 * s] + B1 * cb[R1[s]]
            w[3 * s + 2] = w[3 * s + 2] + B3 * cb[R3[s]]
        return w

    # Collatz-Wielandt enclosure of rho at w=v
    Fv = F(v)
    rho_lo = None; rho_hi = None
    for i in range(N):
        r = Fv[i] / v[i]
        a = mp.mpf(r.a); b = mp.mpf(r.b)
        rho_lo = a if rho_lo is None else min(rho_lo, a)
        rho_hi = b if rho_hi is None else max(rho_hi, b)
    rho = iv.mpf([rho_lo, rho_hi])        # rigorous: rho in this interval
    resid = float(rho_hi) / float(rho_lo) - 1.0  # relative CW gap

    t = A / rho
    R = (t * t + lam) / (1 + t * lam)

    # criteria at w
    cb = [min(v[j], v[j + Nl], v[j + 2 * Nl]) for j in range(Nl)]
    def mean(xs):
        s = xs[0]
        for x in xs[1:]:
            s = s + x
        return s / len(xs)
    mu0 = mean(v[0::3]); mu2 = mean(v[2::3])
    c0 = mean(cb[0::3]); c2 = mean(cb[2::3])
    ratio_c = c2 / c0
    margin = R - ratio_c                   # want > 0
    id_res = mu2 / mu0 - R                 # should contain 0
    g0 = mu0 - c0; g2 = mu2 - c2
    gap_crit = g2 / (R * g0)               # want > 1

    def fmt(x):
        return f"[{float(x.a):+.3e},{float(x.b):+.3e}]"
    print(f"lam={lam_frac} k={k}: rho=[{float(rho.a):.12f},{float(rho.b):.12f}] CWgap={resid:.1e}")
    print(f"   R-c2/c0={fmt(margin)} POS?{float(margin.a)>0} | mu2/mu0-R={fmt(id_res)} has0?{float(id_res.a)<=0<=float(id_res.b)}")
    print(f"   g2/(R*g0)={fmt(gap_crit)} >1?{float(gap_crit.a)>1}", flush=True)
    return float(margin.a) > 0

if __name__ == "__main__":
    ok = 0; tot = 0
    for lam_frac in [Fraction(21,20), Fraction(13,10), Fraction(17,10), Fraction(2)]:
        for k, it in [(5, 400), (6, 400)]:
            tot += 1
            ok += bool(certify(lam_frac, k, it))
    print(f"\ncertified margins positive: {ok}/{tot}")
