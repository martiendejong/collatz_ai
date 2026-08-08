"""
Kernel form of the problem (Obs 528 follow-up):
1. Build the explicit dv/dw kernel via the Neumann sum; verify row-sum identities
   (they must reproduce rho_lin / mass-identity structure in one line).
2. Frozen cb-Jacobian J(j,:) = K_full(j + e*(j)*Nl, :); measure its subleading
   eigenvalue |lam2| (certified-iteration constant; Euler gives lam1 = 1).
"""
import numpy as np
from scipy.sparse import lil_matrix, csr_matrix
from scipy.sparse.linalg import eigs
from math import log2
ALPHA = log2(3.0)

for lam, k in [(1.05, 8), (1.70, 8)]:
    A = lam**-2; B1 = lam**(ALPHA-2); B3 = lam**(ALPHA-1)
    N = 3**(k-1); Nl = N//3
    i = np.arange(N, dtype=np.int64)
    T4 = (4*i+2) % N
    s_arr, r_arr = np.divmod(i, 3)
    m0 = r_arr==0; m2 = r_arr==2
    R1 = (4*s_arr) % Nl; R3 = (2*s_arr+1) % Nl
    v = np.ones(N); rho = 1.0
    for _ in range(2000):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        w = A*v[T4]; w[m2] += B3*cb[R3[m2]]; w[m0] += B1*cb[R1[m0]]
        rho = float(w.max()); w /= rho; v = w
    t = A/rho
    cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])

    # build dv/dw: v(i) = (1/rho) sum_n t^n f(T4^n i); f(i') = B1*w[R1[i']] (class0),
    # B3*w[R3[i']] (class2), 0 (class1)
    nmax = int(np.ceil(np.log(1e-16)/np.log(t)))
    Kf = lil_matrix((N, Nl))
    idx = i.copy(); coef = 1.0/rho
    rowsum_check = np.zeros(N)
    for n in range(nmax):
        cls = (idx % 3)
        s_of = idx // 3
        sel0 = cls == 0; sel2 = cls == 2
        # accumulate coefficients
        rows0 = np.where(sel0)[0]; cols0 = R1[idx[sel0]]
        rows2 = np.where(sel2)[0]; cols2 = R3[idx[sel2]]
        for rr, cc, bb in [(rows0, cols0, B1), (rows2, cols2, B3)]:
            np.add.at(rowsum_check, rr, coef*bb)
        Kf[rows0, cols0] = Kf[rows0, cols0].toarray() + coef*B1 if False else Kf[rows0, cols0]
        # lil fancy add is slow; do dict accumulation instead
        idx = T4[idx]; coef *= t
    # rebuild with COO accumulation (faster correct path)
    from collections import defaultdict
    data = defaultdict(float)
    idx = i.copy(); coef = 1.0/rho
    for n in range(nmax):
        cls = idx % 3
        sel0 = np.where(cls == 0)[0]; sel2 = np.where(cls == 2)[0]
        for rr in sel0:
            data[(rr, R1[idx[rr]])] += coef*B1
        for rr in sel2:
            data[(rr, R3[idx[rr]])] += coef*B3
        idx = T4[idx]; coef *= t
    rows = np.array([p[0] for p in data]); cols = np.array([p[1] for p in data])
    vals = np.array(list(data.values()))
    Kfull = csr_matrix((vals, (rows, cols)), shape=(N, Nl))

    # verify: v = Kfull @ cb (exact reconstruction)
    err = np.abs(Kfull @ cb - v).max()/v.mean()
    # row sums: exact prediction sum over the 3-periodic pattern:
    # row i: sum_n t^n B_{pattern}/rho where pattern depends on class(T4^n i), period 3.
    # class sequence from class(i): c -> 2c+... follows 0->2->1->0; over one period the
    # multiset {B1, B3, 0} always appears once each => rowsum = (B1+B3)(t^{a_i}+...)/..
    # simplest exact check: rowsum = (t^{d0(i)}*Bx + t^{d1(i)}*By)/(rho(1-t^3)) — verify numerically
    rs = np.asarray(Kfull.sum(axis=1)).ravel()
    # per class the rowsum should be constant (only phase differs):
    for r in range(3):
        vals_r = rs[r_arr == r]
        print(f"  lam={lam} class {r}: rowsum mean {vals_r.mean():.8f} sd {vals_r.std():.2e}")
    # predicted from the phase structure: class0 rows start with feed B1 at n=0:
    pred0 = (B1 + B3*t)/(rho*(1-t**3))
    pred2 = (B3 + B1*t**2)/(rho*(1-t**3))
    pred1 = (B1*t + B3*t**2)/(rho*(1-t**3))
    print(f"  voorspeld: klasse0 {pred0:.8f} klasse1 {pred1:.8f} klasse2 {pred2:.8f}")
    print(f"  reconstructie v = K@cb: max rel err {err:.2e}")

    # frozen cb-Jacobian
    trip = np.stack([v[:Nl], v[Nl:2*Nl], v[2*Nl:]], axis=1)
    pi = trip.argmin(axis=1)
    selrows = np.arange(Nl) + pi*Nl
    J = Kfull[selrows, :]
    # Euler: J w* = w*? check
    e1 = np.abs(J @ cb - cb).max()/cb.mean()
    ev = eigs(J, k=4, which='LM', return_eigenvectors=False, maxiter=5000)
    ev = sorted(np.abs(ev))[::-1]
    print(f"  Euler-check |J w* - w*| = {e1:.2e} | |eig| top4: " + " ".join(f"{x:.4f}" for x in ev))
