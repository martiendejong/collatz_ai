"""
Exact fibered (digit-1) mean system:
  mu0[d] = t*mu2[d] + (B1/rho)*c[d]
  mu1[d] = t*mu0[sigma(d)],  sigma(d) = (d+2) mod 3
  mu2[d] = t*mu1[d] + (B3/rho)*c[tau(d)],  tau(d) = (2d+1) mod 3
Cyclic solve: mu0[d] = [f(d) + t^3 f(sigma d) + t^6 f(sigma^2 d)]/(1 - t^9),
  f(d) = t*(B3/rho)*c[tau(d)] + (B1/rho)*c[d].
Validate against measured sub-class means from cached eigenvectors.
"""
import numpy as np
from math import log2

ALPHA = log2(3.0)
CACHE = "E:/projects/collatz/research/cache"

for lam in [1.05, 1.70, 1.90]:
    k = 13
    v = np.load(f"{CACHE}/v_lam{lam:.2f}_k{k}.npy")
    rho = None
    try:
        rho = float(open(f"{CACHE}/rho_lam{lam:.2f}_k{k}.txt").read())
    except Exception:
        # recompute rho quickly via one application
        A_ = lam**-2; B1_ = lam**(ALPHA-2); B3_ = lam**(ALPHA-1)
        N = v.size; Nl = N//3
        i = np.arange(N, dtype=np.int64)
        T4 = (4*i+2) % N
        s_arr, r_arr = np.divmod(i, 3)
        m0 = r_arr==0; m2 = r_arr==2
        R1 = (4*s_arr) % Nl; R3 = (2*s_arr+1) % Nl
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        w = A_*v[T4]; w[m2] += B3_*cb[R3[m2]]; w[m0] += B1_*cb[R1[m0]]
        rho = float((w/v).max())
    A = lam**-2; B1 = lam**(ALPHA-2); B3 = lam**(ALPHA-1)
    t = A/rho
    N = v.size; Nl = N//3
    cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])

    # measured c[d] and sub-class means mu_r[d] (d = sub-class = s mod 3)
    c = np.array([cb[d::3].mean() for d in range(3)])
    s = np.arange(Nl)
    mu_meas = np.zeros((3, 3))
    for r in range(3):
        vr = v[r::3]
        for d in range(3):
            mu_meas[r, d] = vr[s % 3 == d].mean()

    sigma = lambda d: (d+2) % 3
    tau = lambda d: (2*d+1) % 3
    f = np.array([t*(B3/rho)*c[tau(d)] + (B1/rho)*c[d] for d in range(3)])
    mu0 = np.array([(f[d] + t**3*f[sigma(d)] + t**6*f[sigma(sigma(d))])/(1-t**9) for d in range(3)])
    mu1 = np.array([t*mu0[sigma(d)] for d in range(3)])
    mu2 = np.array([t*mu1[d] + (B3/rho)*c[tau(d)] for d in range(3)])
    pred = np.stack([mu0, mu1, mu2])

    err = np.abs(pred - mu_meas).max()/mu_meas.mean()
    print(f"lam={lam}: fibered 9-system max rel err = {err:.2e}")
    # c1 = t*c0 sub-check
    print(f"   c = {c.round(6)}  (c1/(t*c0) = {c[1]/(t*c[0]):.8f}, exact 1?)")
    # digit-1 profile from the system vs measured (overall field)
    prof_pred = pred.mean(axis=0); prof_pred -= prof_pred.mean()
    F = v  # linear-space profile
    i = np.arange(N, dtype=np.int64)
    d1 = (i//3) % 3
    prof_meas = np.array([v[d1==d].mean() for d in range(3)]); prof_meas -= prof_meas.mean()
    cos = prof_pred@prof_meas/np.linalg.norm(prof_pred)/np.linalg.norm(prof_meas)
    print(f"   digit-1 profile: pred {prof_pred.round(6)} meas {prof_meas.round(6)} cos={cos:+.6f}")
