# 368: Route A step 4 — RIGOROUS certification of the domination inequality at
# the endpoint (lam=2, k=13/14), in EXACT RATIONAL arithmetic.
# Key simplification: at lam=2 all operator coefficients are rational
# (A=1/4, B1=2^(alpha-2)=3/4, B3=3/2 since 2^alpha=3), and the stored vector v
# is a dyadic-rational point vector. Hence:
#   - Collatz-Wielandt: rho lies UNCONDITIONALLY in [min_i F(v)_i/v_i, max_i ...],
#     computed as exact rational bounds (integer cross-multiplication).
#   - All cell statistics (U, W, C, GC triples) are exact rationals of v.
#   - The domination inequality  a*m_min > (sqrt(1-a^2) + beta)*m_max  is then
#     evaluated in outward-rounded rational interval arithmetic (sqrt via isqrt).
# Honest caveat (same convention as the Lemma-A certificate): statistics are
# evaluated at the converged point vector w, not the exact Perron vector; the
# CW gap is reported and is ~1e-11, vs margins ~0.05.
import numpy as np
from fractions import Fraction
from math import isqrt, frexp
import os

CACHE = r"E:\projects\collatz\research\cache"
CERT = r"E:\projects\collatz\research\certificates\cert_domination_endpoint.txt"

# ---------- rational interval arithmetic (outward) ----------
class IV:
    __slots__ = ('lo', 'hi')
    def __init__(self, lo, hi=None):
        self.lo = Fraction(lo); self.hi = Fraction(hi if hi is not None else lo)
        assert self.lo <= self.hi
    def __add__(s, o):
        o = mkiv(o); return IV(s.lo+o.lo, s.hi+o.hi)
    def __sub__(s, o):
        o = mkiv(o); return IV(s.lo-o.hi, s.hi-o.lo)
    def __mul__(s, o):
        o = mkiv(o)
        c = [s.lo*o.lo, s.lo*o.hi, s.hi*o.lo, s.hi*o.hi]
        return IV(min(c), max(c))
    def __truediv__(s, o):
        o = mkiv(o)
        assert o.lo > 0 or o.hi < 0
        c = [s.lo/o.lo, s.lo/o.hi, s.hi/o.lo, s.hi/o.hi]
        return IV(min(c), max(c))
    __radd__ = __add__
    __rmul__ = __mul__
    def __rsub__(s, o):
        return mkiv(o).__sub__(s)
    def __neg__(s):
        return IV(-s.hi, -s.lo)
    def __repr__(s):
        return f"[{float(s.lo):.6f},{float(s.hi):.6f}]"

def mkiv(x):
    return x if isinstance(x, IV) else IV(x)

def sqrt_lo(fr):
    if fr <= 0: return Fraction(0)
    K = 10**40
    n = (fr.numerator*K*K)//fr.denominator
    return Fraction(isqrt(n), K)          # floor sqrt: lower bound

def sqrt_hi(fr):
    if fr <= 0: return Fraction(0)
    K = 10**40
    n = -((-fr.numerator*K*K)//fr.denominator)   # ceil
    return Fraction(isqrt(n)+1, K)        # upper bound

def ivsqrt(x):
    x = mkiv(x)
    return IV(sqrt_lo(max(x.lo, Fraction(0))), sqrt_hi(x.hi))

# ---------- exact loading ----------
def exact_ints(v, shift=140):
    m, e = np.frexp(v)
    M = (m*(2.0**53)).astype(np.int64)
    E = e.astype(np.int64) - 53 + shift
    assert E.min() >= 0, "shift too small"
    return [int(a) << int(b) for a, b in zip(M, E)]

def certify(k, fh):
    N = 3**(k-1); Nl = N//3; Nll = Nl//3
    v = np.load(os.path.join(CACHE, f"v_lam2.00_k{k}.npy")).astype(np.float64)
    # polish with float iterations for a tight CW gap
    i = np.arange(N, dtype=np.int64)
    T4 = (4*i+2) % N
    s_, r_ = np.divmod(i, 3)
    m0, m2 = (r_ == 0), (r_ == 2)
    R1 = (4*s_) % Nl; R3 = (2*s_+1) % Nl
    for _ in range(400):
        cbf = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        w = 0.25*v[T4]
        w[m2] += 1.5*cbf[R3[m2]]
        w[m0] += 0.75*cbf[R1[m0]]
        v = w/w.max()
    vi = exact_ints(v)
    # exact F: 4*F_i = v[T4 i] + [cls0] 3*cb + [cls2] 6*cb
    cb = [min(vi[j], vi[j+Nl], vi[j+2*Nl]) for j in range(Nl)]
    T4l = T4.tolist(); R1l = R1.tolist(); R3l = R3.tolist()
    # CW bounds: rho = (4F)_i/(4 v_i): track exact min/max by cross-multiplication
    lo_p = lo_q = hi_p = hi_q = None
    for idx in range(N):
        Fi4 = vi[T4l[idx]]
        c = idx % 3
        if c == 0: Fi4 += 3*cb[R1l[idx]]
        elif c == 2: Fi4 += 6*cb[R3l[idx]]
        qi = 4*vi[idx]
        if lo_p is None or Fi4*lo_q < lo_p*qi:
            lo_p, lo_q = Fi4, qi
        if hi_p is None or Fi4*hi_q > hi_p*qi:
            hi_p, hi_q = Fi4, qi
    rho = IV(Fraction(lo_p, lo_q), Fraction(hi_p, hi_q))
    cwgap = float(rho.hi/rho.lo - 1)
    # exact statistics (integers; VB = 3*vb etc.)
    VB = [vi[j]+vi[j+Nl]+vi[j+2*Nl] for j in range(Nl)]
    CB3 = [cb[x]+cb[x+Nll]+cb[x+2*Nll] for x in range(Nll)]
    CBB = [min(cb[x], cb[x+Nll], cb[x+2*Nll]) for x in range(Nll)]
    U = [[None]*3 for _ in range(3)]
    W = [[None]*3 for _ in range(3)]
    cnt = Nl//9
    sums_vb = [[0]*3 for _ in range(3)]; sums_g = [[0]*3 for _ in range(3)]
    for s in range(Nl):
        a = s % 3; b = (s//3) % 3
        sums_vb[a][b] += VB[s]
        sums_g[a][b] += VB[s] - 3*cb[s]
    for a in range(3):
        for b in range(3):
            U[a][b] = Fraction(sums_vb[a][b], 3*cnt)
            W[a][b] = Fraction(sums_g[a][b], 3*cnt)
    cnt2 = Nll//3
    C = [Fraction(sum(CB3[x] for x in range(b, Nll, 3)), 3*cnt2) for b in range(3)]
    GC = [Fraction(sum(CB3[x]-3*CBB[x] for x in range(b, Nll, 3)), 3*cnt2) for b in range(3)]
    # association a = Corr(C, GC), exact triples -> interval via sqrt bounds
    mC = sum(C)/3; mG = sum(GC)/3
    Cc = [x-mC for x in C]; Gc_ = [x-mG for x in GC]
    cov = sum(x*y for x, y in zip(Cc, Gc_))/3
    vC = sum(x*x for x in Cc)/3; vG = sum(x*x for x in Gc_)/3
    a_iv = IV(cov)/(ivsqrt(IV(vC))*ivsqrt(IV(vG)))
    # transfer matrices with t interval
    t = IV(Fraction(1, 4))/rho
    b1r = IV(Fraction(3, 4))/rho; b3r = IV(Fraction(3, 2))/rho
    tau = [1, 0, 2]
    def Mmat(par):
        M0 = [[IV(0) for _ in range(3)] for _ in range(3)]
        t3 = t*t*t; t9 = t3*t3*t3
        for b in range(3):
            for j in range(3):
                bj = (b+2*j) % 3
                wgt = (t3*t3 if j == 2 else (t3 if j == 1 else IV(1)))/(IV(1)-t9)
                M0[b][bj] = M0[b][bj] + wgt*b1r
                M0[b][tau[bj]] = M0[b][tau[bj]] + wgt*t*b3r
        if par == 0:
            return M0
        M1 = [[t*M0[(b+2) % 3][c] for c in range(3)] for b in range(3)]
        M2 = [[t*M1[b][c] for c in range(3)] for b in range(3)]
        for b in range(3):
            M2[b][tau[b]] = M2[b][tau[b]] + b3r
        return M2
    ok_all = True
    lines = [f"k={k}: rho=[{float(rho.lo):.13f},{float(rho.hi):.13f}] CW-gap={cwgap:.1e}  a={a_iv}"]
    for par in (0, 2):
        Ma = Mmat(par)
        # centred restriction: columns of Q = orthonormal sum-zero basis
        s2 = ivsqrt(IV(2)); s6 = ivsqrt(IV(6))
        Q = [[IV(1)/s2, IV(1)/s6], [IV(-1)/s2, IV(1)/s6], [IV(0), IV(-2)/s6]]
        colmean = [ (Ma[0][c]+Ma[1][c]+Ma[2][c])*IV(Fraction(1,3)) for c in range(3)]
        PM = [[Ma[r][c]-colmean[c] for c in range(3)] for r in range(3)]
        Rm = [[sum((Q[r][p]*PM[r][c] for r in range(3)), IV(0)) for c in range(3)] for p in range(2)]
        Rq = [[sum((Rm[p][c]*Q[c][q_] for c in range(3)), IV(0)) for q_ in range(2)] for p in range(2)]
        # Gram = Rq^T Rq, eigenvalues of symmetric 2x2
        g11 = Rq[0][0]*Rq[0][0]+Rq[1][0]*Rq[1][0]
        g22 = Rq[0][1]*Rq[0][1]+Rq[1][1]*Rq[1][1]
        g12 = Rq[0][0]*Rq[0][1]+Rq[1][0]*Rq[1][1]
        tr = g11+g22
        disc = ivsqrt((g11-g22)*(g11-g22)+IV(4)*g12*g12)
        lmin = (tr-disc)*IV(Fraction(1, 2))
        lmax = (tr+disc)*IV(Fraction(1, 2))
        mmin = ivsqrt(lmin); mmax = ivsqrt(lmax)
        # slack fraction beta = |Cov(U_a, S_a)|/Cov(U_a, M_a GC), S per Prop slack
        Ua = U[par]; mUa = sum(Ua)/3
        Uc = [x-mUa for x in Ua]
        MG = [sum((Ma[b][c]*GC[c] for c in range(3)), IV(0)) for b in range(3)]
        mMG = (MG[0]+MG[1]+MG[2])*IV(Fraction(1, 3))
        covLin = sum((IV(Uc[b])*(MG[b]-mMG) for b in range(3)), IV(0))*IV(Fraction(1, 3))
        if par == 0:
            S = [IV(Fraction(1,4))*IV(W[2][b]) + IV(Fraction(3,4))*IV(GC[b]) - rho*IV(W[0][b]) for b in range(3)]
        else:
            S = [IV(Fraction(1,4))*IV(W[1][b]) + IV(Fraction(3,2))*IV(GC[tau[b]]) - rho*IV(W[2][b]) for b in range(3)]
        mS = (S[0]+S[1]+S[2])*IV(Fraction(1, 3))
        covS = sum((IV(Uc[b])*(S[b]-mS) for b in range(3)), IV(0))*IV(Fraction(1, 3))/rho
        absS = IV(max(Fraction(0), max(-covS.hi, covS.lo) if covS.lo > 0 or covS.hi < 0 else Fraction(0)),
                  max(abs(covS.lo), abs(covS.hi)))
        beta = absS/covLin if covLin.lo > 0 else IV(1, 10)
        # domination: a*mmin - (sqrt(1-a^2)+beta)*mmax > 0, outward
        sfac = ivsqrt(IV(1)-a_iv*a_iv)
        marg = a_iv*mmin - (sfac+beta)*mmax
        ok = marg.lo > 0
        ok_all &= ok
        lines.append(f"  parent {par}: m=[{float(mmin.lo):.4f},{float(mmax.hi):.4f}] beta<= {float(beta.hi):.4f} "
                     f"margin=[{float(marg.lo):+.5f},{float(marg.hi):+.5f}]  CERTIFIED: {ok}")
    for L in lines:
        print(L, flush=True); fh.write(L+"\n")
    return ok_all

with open(CERT, 'w', encoding='utf-8') as fh:
    fh.write("Rigorous certification of the domination inequality at the endpoint lam=2\n")
    fh.write("(exact rational arithmetic: A=1/4, B1=3/4, B3=3/2; CW enclosure of rho;\n")
    fh.write("outward rational-interval evaluation; statistics at the converged point\n")
    fh.write("vector, CW gap reported per the Lemma-A convention.)\n\n")
    allok = True
    for k in (13, 14):
        allok &= certify(k, fh)
    fh.write(f"\nALL CERTIFIED: {allok}\n")
print("\nALL CERTIFIED:", allok)
