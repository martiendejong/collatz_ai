"""Obs 625: P2 · search for a Lyapunov metric for the K-L linearization.

Goal: find diagonal Q = diag(C^s) such that the linearization M at the
true Perron fixed point, restricted to the fine-mode space
S = {x : u.x = 0} (M-invariant), has Q-operator norm < 1 · and check
whether the optimal s is stable across k. A k-stable metric with norm
uniformly < 1 converts the Open Lemma into an induction step over the
exact lattice identity (Thm 16).

Also reports the spectral radius of M on S, which should reproduce the
measured attenuation kappa ~ 0.84 (validation of the whole framing).

Method: (1) iterate the nonlinear operator T(c) = W0 c(i4) + Wb min-triple
from the certificate to the true eigenpair; (2) freeze the argmin
pattern -> sparse M; (3) left vector u by transpose iteration;
(4) for each s: Q-orthogonal projection onto S, power iteration on
P M*_Q M P for the norm; plain deflated iteration for the radius.
"""
import numpy as np

def build(K, lam, path, t_iters=300):
    N = 3 ** (K - 1); M3 = 3 ** K; Mc = 3 ** (K - 1)
    beta = np.log2(3); W0 = lam ** -2; W2 = lam ** (beta - 2); W8 = lam ** (beta - 1)
    C = np.asarray(np.load(path, mmap_mode='r')[:]).astype(np.float64)
    idx = np.arange(N, dtype=np.int64); m = 3 * idx + 2
    i4 = (((4 * m) % M3) - 2) // 3; mod9 = m % 9
    J = {}; Wb = np.zeros(N)
    for mask, mul, r, W in ((mod9 == 2, 4, 2, W2), (mod9 == 8, 2, 1, W8)):
        mm = m[mask]; t = ((mul * mm - r) // 3) % Mc
        J[W] = (mask, np.stack([(t - 2) // 3, ((t + Mc) - 2) // 3, ((t + 2 * Mc) - 2) // 3]))
        Wb[mask] = W
    c = C / C.sum()
    for it in range(t_iters):
        Tc = W0 * c[i4]
        for W, (mask, j) in J.items():
            Tc[mask] += W * np.min(c[j], axis=0)
        lam_est = Tc.sum() / c.sum()
        Tc /= Tc.sum()
        if np.abs(Tc - c).sum() < 1e-13:
            c = Tc; break
        c = Tc
    tri = np.full(N, -1, dtype=np.int64)
    for W, (mask, j) in J.items():
        o = np.argmin(c[j], axis=0)
        tri[mask] = j[o, np.arange(j.shape[1])]
    hm = tri >= 0
    def Mx(x):
        y = W0 * x[i4]
        y[hm] += Wb[hm] * x[tri[hm]]
        return y
    def MTx(y):
        z = np.bincount(i4, weights=W0 * y, minlength=N)
        z += np.bincount(tri[hm], weights=Wb[hm] * y[hm], minlength=N)
        return z
    u = np.full(N, 1.0 / N)
    for it in range(300):
        nu = MTx(u); mu = nu.sum(); nu /= nu.sum()
        if np.abs(nu - u).sum() < 1e-13:
            u = nu; break
        u = nu
    print(f"k={K}: nonlinear eigen-iteration lambda_norm = {lam_est:.6f} (want ~1); "
          f"argmin rows {hm.sum():,}")
    return N, c, u, Mx, MTx

def fine_radius(N, u, Mx, iters=400, seed=3):
    rng = np.random.default_rng(seed)
    x = rng.standard_normal(N)
    x -= u * (u @ x) / (u @ u)
    r_prev = 0.0
    for it in range(iters):
        y = Mx(x)
        y -= u * (u @ y) / (u @ u)
        r = np.linalg.norm(y) / np.linalg.norm(x)
        x = y / np.linalg.norm(y)
        if it > 50 and abs(r - r_prev) < 1e-8:
            break
        r_prev = r
    return r

def fine_qnorm(N, c, u, Mx, MTx, s, iters=250, seed=5):
    q = np.power(c / c.mean(), s)
    qinv_u = u / q
    denom = u @ qinv_u
    def P(x):
        return x - qinv_u * (u @ x) / denom
    rng = np.random.default_rng(seed)
    z = P(rng.standard_normal(N))
    nrm = None
    for it in range(iters):
        w = Mx(P(z))
        w2 = MTx(q * w) / q
        z2 = P(w2)
        nrm_new = np.sqrt(max((z2 * q * z2).sum(), 0) ) / np.sqrt((z * q * z).sum())
        z = z2 / np.linalg.norm(z2)
        if nrm is not None and abs(nrm_new - nrm) < 1e-7 and it > 60:
            nrm = nrm_new; break
        nrm = nrm_new
    return np.sqrt(nrm)

if __name__ == "__main__":
    for K, lam, path in ((13, 1.818, r"E:\projects\collatz\research\certificates\cert_k13.npy"),
                         (15, 1.841, r"E:\projects\collatz\research\certificates\cert_k15.npy")):
        N, c, u, Mx, MTx = build(K, lam, path)
        rad = fine_radius(N, u, Mx)
        print(f"  spectral radius on fine space: {rad:.4f}  (measured kappa ~ 0.839)")
        best = (None, np.inf)
        for s in (-1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0):
            nrm = fine_qnorm(N, c, u, Mx, MTx, s)
            tag = " <-- contraction metric" if nrm < 1 else ""
            print(f"  s={s:+.1f}: fine-space Q-norm = {nrm:.4f}{tag}")
            if nrm < best[1]:
                best = (s, nrm)
        print(f"  best: s={best[0]:+.1f}, norm {best[1]:.4f}\n")
