# 366: Route A step 2 — the transfer algebra of the q=1 rung (follow-up Obs 552).
# Exact claims to verify at machine precision, then quantify the rest:
# (T1) cell-mean 9-system for U (from the exact tower identity + affine cell
#      bijections):  rho*U[0,b] = A*U[2,b] + B1*C[b]
#                    rho*U[1,b] = A*U[0,(b+2)%3]
#                    rho*U[2,b] = A*U[1,b] + B3*C[tau(b)],  tau = (0 1)(2)
#      hence U = M*C with an explicit POSITIVE matrix (denominator 1-t^9).
# (T2) NEW gap-copy identity: on class-1 fibers both vb and cb copy pointwise
#      (v is a t-copy on class 1, and all lifts stay class 1), so
#      G(s) = t*G(sigma4 s) for s = 1 mod 3, hence W[1,b] = t*W[0,(b+2)%3]
#      and Cov(parent 1) = t^2 * Cov(parent 0): parent 1 needs no proof.
# (T3) min-superadditivity at fiber level gives one-sided linear bounds:
#      rho*W[0,b] <= A*W[2,b] + B1*GC[b]   (slack s1[b] >= 0)
#      rho*W[2,b] <= A*W[1,b] + B3*GC[tau(b)] (slack s3[b] >= 0)
# (Q) quantify: solve the slackless system What = M*GC; compare Cov(U, What)
#     with the true Cov(U, W) per parent: the deficit is the slack-carried part
#     (Lemma-beta territory). If the linear part dominates, the sign of b1
#     rides on the input association (measured 0.98, Obs 552).
import numpy as np
import os
from math import log2

ALPHA = log2(3.0)
CACHE = r"E:\projects\collatz\research\cache"

def load_or_make(lam, k, iters=1200):
    N = 3**(k-1)
    fn = os.path.join(CACHE, f"v_lam{lam:.2f}_k{k}.npy")
    i = np.arange(N, dtype=np.int64)
    T4 = (4*i+2) % N
    s_, r_ = np.divmod(i, 3)
    Nl = N//3
    m0, m2 = (r_ == 0), (r_ == 2)
    R1 = (4*s_) % Nl; R3 = (2*s_+1) % Nl
    A = lam**-2; B1 = lam**(ALPHA-2); B3 = lam**(ALPHA-1)
    v = np.load(fn).astype(np.float64) if os.path.exists(fn) else np.ones(N)
    it = 300 if os.path.exists(fn) else iters
    rho = 1.0
    for _ in range(it):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        w = A*v[T4]
        w[m2] += B3*cb[R3[m2]]
        w[m0] += B1*cb[R1[m0]]
        rho = float(w.max()); v = w/rho
    return v, rho, A, B1, B3

def M_matrix(t, b1r, b3r):
    # solve the 9-system symbolically: U0[b] = t^3 U0[b+2] + t*b3r*C[tau b] + b1r*C[b]
    # over the 3-cycle b -> b+2 -> b+1 -> b. Returns 3x3 matrices M0, M1, M2 with
    # U_a = M_a @ C.
    tau = [1, 0, 2]
    M0 = np.zeros((3, 3))
    for b in range(3):
        # iterate j = 0,1,2 with b_j = (b + 2j) % 3, weight t^{3j}, divided by (1-t^9)
        for j in range(3):
            bj = (b + 2*j) % 3
            w = t**(3*j)/(1 - t**9)
            M0[b, bj] += w*b1r
            M0[b, tau[bj]] += w*t*b3r
    M1 = np.zeros((3, 3))
    M2 = np.zeros((3, 3))
    for b in range(3):
        M1[b] = t*M0[(b+2) % 3]                      # U1[b] = t U0[b+2]
    for b in range(3):
        M2[b] = t*M1[b]                              # U2[b] = t U1[b] + b3r C[tau b]
        M2[b, tau[b]] += b3r
    return M0, M1, M2

def analyse(lam, k):
    v, rho, A, B1, B3 = load_or_make(lam, k)
    N = v.size; Nl = N//3; Nll = Nl//3
    s = np.arange(Nl, dtype=np.int64)
    vb = (v[:Nl] + v[Nl:2*Nl] + v[2*Nl:])/3.0
    cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
    G = vb - cb
    cb3 = (cb[:Nll] + cb[Nll:2*Nll] + cb[2*Nll:])/3.0
    cbb = np.minimum(np.minimum(cb[:Nll], cb[Nll:2*Nll]), cb[2*Nll:])
    Gc = cb3 - cbb
    d0 = s % 3; d1 = (s // 3) % 3
    U = np.zeros((3, 3)); W = np.zeros((3, 3))
    for a in range(3):
        for b in range(3):
            m = (d0 == a) & (d1 == b)
            U[a, b] = vb[m].mean(); W[a, b] = G[m].mean()
    x2 = np.arange(Nll)
    C = np.array([cb3[x2 % 3 == b].mean() for b in range(3)])
    GC = np.array([Gc[x2 % 3 == b].mean() for b in range(3)])
    t = A/rho; b1r = B1/rho; b3r = B3/rho
    tau = [1, 0, 2]
    # T1: 9-system residuals + closed form
    r1 = max(abs(rho*U[0, b] - A*U[2, b] - B1*C[b]) for b in range(3))
    r2 = max(abs(rho*U[1, b] - A*U[0, (b+2) % 3]) for b in range(3))
    r3 = max(abs(rho*U[2, b] - A*U[1, b] - B3*C[tau[b]]) for b in range(3))
    M0, M1, M2 = M_matrix(t, b1r, b3r)
    rM = max(np.abs(M0@C - U[0]).max(), np.abs(M1@C - U[1]).max(), np.abs(M2@C - U[2]).max())
    # T2: pointwise gap copy on class-1 fibers
    sig4 = (4*s+2) % Nl
    m1f = (d0 == 1)
    gcopy = float(np.abs(G[m1f] - t*G[sig4[m1f]]).max()/(np.abs(G).max()+1e-300))
    wcopy = max(abs(W[1, b] - t*W[0, (b+2) % 3]) for b in range(3))
    covs = [float(np.cov(U[a], W[a], bias=True)[0, 1]) for a in range(3)]
    cov_ratio = covs[1]/(t*t*covs[0])
    # T3: slack signs
    s1 = [A*W[2, b] + B1*GC[b] - rho*W[0, b] for b in range(3)]
    s3 = [A*W[1, b] + B3*GC[tau[b]] - rho*W[2, b] for b in range(3)]
    # Q: slackless covariance vs true
    What0 = M0@GC; What2 = M2@GC
    cov_lin0 = float(np.cov(U[0], What0, bias=True)[0, 1])
    cov_lin2 = float(np.cov(U[2], What2, bias=True)[0, 1])
    return dict(r_sys=max(r1, r2, r3), r_M=rM, gcopy=gcopy, wcopy=wcopy,
                cov_ratio=cov_ratio, s1=min(s1), s3=min(s3),
                lin_frac0=cov_lin0/covs[0] if covs[0] else np.nan,
                lin_frac2=cov_lin2/covs[2] if covs[2] else np.nan,
                scale=abs(U).max())

print(f"{'lam':>5} {'k':>3} {'9sys-res':>9} {'U=MC-res':>9} {'G-copy':>8} {'W-copy':>8} "
      f"{'cov1/t2cov0':>11} {'min s1':>9} {'min s3':>9} {'lin0':>6} {'lin2':>6}")
allok = True
for lam in [1.05, 1.30, 1.70, 2.00]:
    for k in [10, 13, 14]:
        d = analyse(lam, k)
        sc = d['scale']
        print(f"{lam:>5} {k:>3} {d['r_sys']/sc:>9.1e} {d['r_M']/sc:>9.1e} {d['gcopy']:>8.1e} "
              f"{d['wcopy']/sc:>8.1e} {d['cov_ratio']:>11.4f} {d['s1']/sc:>9.2e} {d['s3']/sc:>9.2e} "
              f"{d['lin_frac0']:>6.3f} {d['lin_frac2']:>6.3f}", flush=True)
        if d['s1'] < -1e-12*sc or d['s3'] < -1e-12*sc:
            allok = False
print()
print("alle slacks niet-negatief (T3):", allok)
