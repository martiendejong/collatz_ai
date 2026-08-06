"""
Birkhoff step, quantitative skeleton (case lam=21/20, k=7 — smallest margin).
1. Converge w at high precision; measure argmin STABILITY margin:
   min over cb-columns of (second_smallest - smallest)/scale.
   If localization radius << this margin, F == M_pi (linear) near w and v*.
2. Build sparse M_pi; find positivity depth m (M_pi^m > 0 entrywise).
3. Bound projective diameter Delta_m = max_{i,j} d_H(col_i, col_j) of M_pi^m;
   kappa = tanh(Delta_m/4); localization d_H(w, v*) <= m*eps/(1-kappa).
"""
import numpy as np
from math import log2, tanh, log

ALPHA = log2(3.0)
lam = 21/20
k = 7
A = lam**-2; B1 = lam**(ALPHA-2); B3 = lam**(ALPHA-1)
N = 3**(k-1); Nl = N//3
i = np.arange(N, dtype=np.int64)
T4 = (4*i+2) % N
s_arr, r_arr = np.divmod(i, 3)
m0 = r_arr==0; m2 = r_arr==2
R1 = (4*s_arr) % Nl; R3 = (2*s_arr+1) % Nl
v = np.ones(N)
for _ in range(3000):
    cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
    w_ = A*v[T4]
    w_[m2] += B3*cb[R3[m2]]
    w_[m0] += B1*cb[R1[m0]]
    rho = w_.max(); v = w_/rho

# 1. argmin stability margin
trip = np.stack([v[:Nl], v[Nl:2*Nl], v[2*Nl:]], axis=1)
srt = np.sort(trip, axis=1)
stab = ((srt[:,1]-srt[:,0])/srt[:,0])
print(f"argmin stability: min relative margin = {stab.min():.3e} (median {np.median(stab):.3e})")

# 2. frozen linear operator M_pi (sparse, 2 entries/row)
pi = trip.argmin(axis=1)
rows = []; cols = []; vals = []
for ii in range(N):
    rows.append(ii); cols.append(T4[ii]); vals.append(A)
for s in range(Nl):
    j = R1[3*s]; rows.append(3*s); cols.append(j + pi[j]*Nl); vals.append(B1)
    j3 = R3[3*s]; rows.append(3*s+2); cols.append(j3 + pi[j3]*Nl); vals.append(B3)
from scipy.sparse import csr_matrix
M = csr_matrix((vals, (rows, cols)), shape=(N,N))
# fixed point check
resid = np.abs(M@v - rho*v).max()/v.min()
print(f"frozen-linear check: |M w - rho w|_inf / min(w) = {resid:.2e}")

# positivity depth
P = M.copy(); P.data[:] = 1.0
depth = 1
X = P.copy()
while (X.count_nonzero() < N*N) and depth < 60:
    X = (X @ P)
    X.data[:] = 1.0
    depth += 1
print(f"positivity depth m = {depth} (full: {X.count_nonzero()==N*N})")

# 3. projective diameter of M^m: Delta = max_{j,j'} max_i log((c_j[i]/c_j'[i]) * (c_j'[l]/c_j[l]))
m = depth
Mm = M.copy()
for _ in range(m-1):
    Mm = Mm @ M
Mm = np.asarray(Mm.todense(), dtype=float)
if (Mm > 0).all():
    # Delta = max over column pairs of max_i,l log( Mm[i,j]*Mm[l,j'] / (Mm[l,j]*Mm[i,j']) )
    # = max_{j,j'} [ max_i log(Mm[i,j]/Mm[i,j']) - min_i log(...) ]
    L = np.log(Mm)
    # for each pair of columns: range of (L[:,j]-L[:,j'])
    D = 0.0
    for j in range(N):
        diff = L - L[:, j][:, None]   # (i, j') matrix: L[i,j'] - L[i,j]
        rng = diff.max(axis=0) - diff.min(axis=0)
        D = max(D, rng.max())
    kappa = tanh(D/4)
    eps = 1e-30   # CW residual scale from certification
    loc = m*eps/(1-kappa) if kappa < 1 else float('inf')
    print(f"Delta_{m} = {D:.3f}  kappa = {kappa:.6f}  localization d_H(w,v*) <= {loc:.2e}")
    print(f"vs argmin stability {stab.min():.2e} and criterion margin 5.4e-4: "
          f"{'OK — chain closes' if loc < 1e-10 else 'insufficient'}")
