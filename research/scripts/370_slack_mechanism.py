# 370: Route A step 6 — the MECHANISM of the slack fraction beta.
# Questions:
#  (S1) do the slack triples S1[b], S3[b] inherit structure (copy relations,
#       shape alignment with U or GC)?
#  (S2) does beta DECAY with k (367 hinted: 0.184 -> 0.143 at the endpoint),
#       and at which rate? Candidates: sqrt(c) (amplitude rate) or the margin
#       rate. If beta -> 0 geometrically, the domination inequality becomes
#       asymptotically automatic and the whole induction hypothesis reduces to
#       the margin law alone.
#  (S3) second-order hypothesis: slack ~ misalignment cost ~ local fluctuation
#       (second order in the sliver), so slack/level should decay at the
#       amplitude rate; measure the relative slack scale.
import numpy as np
import os
from math import log2

ALPHA = log2(3.0)
CACHE = r"E:\projects\collatz\research\cache"

def load_or_make(lam, k, iters=1400):
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
    it = 400 if os.path.exists(fn) else iters
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
    return M0, M2

def cs(u, w):
    uc = u - u.mean(); wc = w - w.mean()
    n = np.sqrt((uc@uc)*(wc@wc))
    return float(uc@wc/n) if n > 0 else 0.0

rows = {}
print(f"{'lam':>5} {'k':>3} {'S1/S1[0]':>22} {'S3/S3[0]':>22} {'S1~U':>6} {'S1~GC':>6} "
      f"{'beta0':>7} {'beta2':>7} {'relS':>9}")
for lam in [1.05, 1.30, 1.70, 2.00]:
    for k in [10, 12, 13, 14]:
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
        GC = np.array([Gc[x2 % 3 == b_].mean() for b_ in range(3)])
        t = A/rho
        tau = [1, 0, 2]
        S1 = np.array([A*W[2, b_] + B1*GC[b_] - rho*W[0, b_] for b_ in range(3)])
        S3 = np.array([A*W[1, b_] + B3*GC[tau[b_]] - rho*W[2, b_] for b_ in range(3)])
        M0, M2 = M_matrices(t, B1/rho, B3/rho)
        beta = []
        for par, Ma, S in [(0, M0, S1), (2, M2, S3)]:
            Ua = U[par] - U[par].mean()
            MG = Ma@GC; MG = MG - MG.mean()
            covLin = float(np.mean(Ua*MG))
            Sc = S - S.mean()
            covS = float(np.mean(Ua*Sc))/rho
            beta.append(abs(covS)/abs(covLin))
        relS = S1.mean()/ (rho*W[0].mean())
        rows.setdefault(lam, []).append((k, beta[0], beta[1], relS))
        print(f"{lam:>5} {k:>3} {np.array2string(S1/S1[0], precision=3):>22} "
              f"{np.array2string(S3/S3[0], precision=3):>22} {cs(S1, U[0]):>6.3f} {cs(S1, GC):>6.3f} "
              f"{beta[0]:>7.4f} {beta[1]:>7.4f} {relS:>9.5f}", flush=True)
print()
print("verval per dieptestap:")
from math import sqrt
c_est = {1.05: 0.41, 1.30: 0.55, 1.70: 0.70, 2.00: 0.835}
for lam, rr in rows.items():
    rr.sort()
    ks = [r[0] for r in rr]
    b0 = [r[1] for r in rr]; b2 = [r[2] for r in rr]; rl = [r[3] for r in rr]
    rb0 = [(b0[i+1]/b0[i])**(1/(ks[i+1]-ks[i])) for i in range(len(rr)-1)]
    rb2 = [(b2[i+1]/b2[i])**(1/(ks[i+1]-ks[i])) for i in range(len(rr)-1)]
    rrl = [(rl[i+1]/rl[i])**(1/(ks[i+1]-ks[i])) for i in range(len(rr)-1)]
    print(f"  lam={lam}: rate(beta0) {['%.3f' % r for r in rb0]}  rate(beta2) {['%.3f' % r for r in rb2]}  "
          f"rate(relS) {['%.3f' % r for r in rrl]}  sqrt(c)~{sqrt(c_est[lam]):.3f}")
