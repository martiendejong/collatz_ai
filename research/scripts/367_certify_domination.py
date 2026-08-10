# 367: Route A step 3 — make f(beta*) a NUMBER: explicit singular values of the
# transfer matrices M_a on the centred-triple subspace, and per-grid-point check
# of the domination condition of the conditional induction lemma:
#   conservative:  a*m_min > (sqrt(1-a^2) + beta) * m_max
#   exact:         Cov(U_a, W_a) = kappa*Cov(U_a, M_a C) + Cov(U_a, M_a eps) - Cov(U_a, S_a)
# We report both: where the conservative form certifies, and where only the
# exact decomposition does (honest gap for the write-up).
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

def M_matrices(t, b1r, b3r):
    tau = [1, 0, 2]
    M0 = np.zeros((3, 3))
    for b in range(3):
        for j in range(3):
            bj = (b + 2*j) % 3
            w = t**(3*j)/(1 - t**9)
            M0[b, bj] += w*b1r
            M0[b, tau[bj]] += w*t*b3r
    M1 = np.array([t*M0[(b+2) % 3] for b in range(3)])
    M2 = np.array([t*M1[b] for b in range(3)])
    for b in range(3):
        M2[b, tau[b]] += b3r
    return M0, M1, M2

def centred_sv(M):
    # restrict to the 2-dim sum-zero subspace: orthonormal basis Q (3x2)
    Q = np.array([[1, -1, 0], [1, 1, -2]], dtype=float).T
    Q[:, 0] /= np.linalg.norm(Q[:, 0]); Q[:, 1] /= np.linalg.norm(Q[:, 1])
    # centred action: P M Q where P projects output to sum-zero then expresses in Q
    PM = M - M.mean(axis=0, keepdims=True)   # output centring
    R = Q.T @ PM @ Q
    sv = np.linalg.svd(R, compute_uv=False)
    return sv.min(), sv.max()

print(f"{'lam':>5} {'k':>3} {'par':>3} {'m_min':>8} {'m_max':>8} {'cond':>6} {'a':>7} {'beta':>7} "
      f"{'cons.margin':>11} {'OK?':>4} {'exact Cov>0':>11}")
n_cons_ok = n_exact_ok = n_tot = 0
for lam in [1.05, 1.30, 1.70, 2.00]:
    for k in [10, 13, 14]:
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
        for a_ in range(3):
            for b_ in range(3):
                m = (d0 == a_) & (d1 == b_)
                U[a_, b_] = vb[m].mean(); W[a_, b_] = G[m].mean()
        x2 = np.arange(Nll)
        C = np.array([cb3[x2 % 3 == b_].mean() for b_ in range(3)])
        GC = np.array([Gc[x2 % 3 == b_].mean() for b_ in range(3)])
        t = A/rho
        M0, M1, M2 = M_matrices(t, B1/rho, B3/rho)
        Cc = C - C.mean(); GCc = GC - GC.mean()
        a_assoc = float(Cc @ GCc/(np.linalg.norm(Cc)*np.linalg.norm(GCc) + 1e-300))
        tau = [1, 0, 2]
        for par, Ma in [(0, M0), (2, M2)]:
            mmin, mmax = centred_sv(Ma)
            Ua = U[par] - U[par].mean()
            Wa = W[par] - W[par].mean()
            covUW = float(np.mean(Ua*Wa))
            covLin = float(np.mean(Ua*(Ma@GC - (Ma@GC).mean())))
            # slack per cell
            if par == 0:
                S = np.array([A*W[2, b_] + B1*GC[b_] - rho*W[0, b_] for b_ in range(3)])
            else:
                S = np.array([A*W[1, b_] + B3*GC[tau[b_]] - rho*W[2, b_] for b_ in range(3)])
            Sc = S - S.mean()
            covS = float(np.mean(Ua*Sc))/rho
            beta = abs(covS)/(abs(covLin) + 1e-300)
            cons = a_assoc*mmin - (np.sqrt(max(0, 1-a_assoc**2)) + beta)*mmax
            ok = cons > 0
            n_tot += 1
            n_cons_ok += ok
            n_exact_ok += (covUW > 0)
            print(f"{lam:>5} {k:>3} {par:>3} {mmin:>8.4f} {mmax:>8.4f} {mmax/mmin:>6.2f} "
                  f"{a_assoc:>7.4f} {beta:>7.3f} {cons:>11.4f} {'JA' if ok else 'nee':>4} "
                  f"{'+' if covUW > 0 else '-':>11}", flush=True)
print()
print(f"conservatieve dominantie gecertificeerd: {n_cons_ok}/{n_tot}")
print(f"exacte covariantie positief:             {n_exact_ok}/{n_tot}")
