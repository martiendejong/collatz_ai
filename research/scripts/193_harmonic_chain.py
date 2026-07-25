"""
193_harmonic_chain.py
=====================
THE (S)-ARCHITECTURE TEST: the edge equation defines an exact Markov chain
("harmonic chain"): P(m -> 4m) = lam^-2 v(4m)/v(m), P(m -> sigma(m)) =
phi(m); rows sum to the edge growth ~ 1 (fixed-point identity).

Architecture (to be tested):
 (1) rows sum to 1 at the edge                        [exact, +polish residue]
 (2) under the invariant measure pi ~ v*u (u = left Perron vector), the
     feed-time subsequence is stationary => theta_g constant in g
 (3) theta* < 1 quantitatively (two-sided field bounds on v AND u)
 (4) mu-start transient controlled by the density sandwich d(mu)/d(pi)

Measurements at k=13: row sums; Var(log u); corr(log v, log u); density
spread mu/pi; theta_g cascade started from mu vs started from pi.
"""
import numpy as np
from math import log2

ALPHA = log2(3.0)

def make_maps(k):
    N = 3**(k-1)
    i = np.arange(N, dtype=np.int64)
    T4 = (4*i+2) % N
    s, r = np.divmod(i, 3)
    Nl = N//3
    return N, Nl, T4, (r == 0), (4*s) % Nl, (r == 2), (2*s+1) % Nl

def edge(k):
    N, Nl, T4, m1, R1, m3, R3 = make_maps(k)
    lo, hi = 1.5, 1.999
    v = np.ones(N)
    for _ in range(36):
        lam = (lo+hi)/2
        A, B1, B3 = lam**-2, lam**(ALPHA-2), lam**(ALPHA-1)
        w = v.copy()
        for _ in range(60):
            cb = np.minimum(np.minimum(w[:Nl], w[Nl:2*Nl]), w[2*Nl:])
            w2 = A*w[T4]
            w2[m1] += B1*cb[R1[m1]]
            w2[m3] += B3*cb[R3[m3]]
            g = w2.max()
            w = w2/g
        if g >= 1.0:
            lo, v = lam, w
        else:
            hi = lam
    lam = lo
    A, B1, B3 = lam**-2, lam**(ALPHA-2), lam**(ALPHA-1)
    for _ in range(400):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        w2 = A*v[T4]
        w2[m1] += B1*cb[R1[m1]]
        w2[m3] += B3*cb[R3[m3]]
        v = w2/w2.max()
    return lam, v

k = 13
lam, v = edge(k)
N, Nl, T4, m1, R1, m3, R3 = make_maps(k)
A, B1, B3 = lam**-2, lam**(ALPHA-2), lam**(ALPHA-1)

stack = np.stack([v[:Nl], v[Nl:2*Nl], v[2*Nl:]])
amin = stack.argmin(axis=0)
rf = np.zeros(N, dtype=np.int64)
rf[m1] = R1[m1]
rf[m3] = R3[m3]
has = m1 | m3
sigma = np.full(N, -1, dtype=np.int64)
sigma[has] = rf[has] + amin[rf[has]]*Nl
bcoef = np.zeros(N)
bcoef[m1] = B1
bcoef[m3] = B3

p_back = A*v[T4]/v
p_feed = np.zeros(N)
p_feed[has] = bcoef[has]*v[sigma[has]]/v[has]
rowsum = p_back + p_feed
print(f"k={k} lam={lam:.4f}")
print(f"(1) harmonic row sums: mean {rowsum.mean():.6f} "
      f"min {rowsum.min():.6f} max {rowsum.max():.6f}")

# left Perron vector of the frozen system
invT4 = np.empty(N, dtype=np.int64)
invT4[T4] = np.arange(N)
u = np.ones(N)
for _ in range(400):
    un = A*u[invT4]
    np.add.at(un, sigma[has], bcoef[has]*u[has])
    u = un/un.max()
pi = v*u
pi /= pi.sum()
mu = v/v.sum()
lu = np.log2(u)
lv = np.log2(v)
print(f"(3) Var(log2 u) = {lu.var():.4f}   Var(log2 v) = {lv.var():.4f}   "
      f"corr(log v, log u) = {np.corrcoef(lv, lu)[0,1]:+.4f}")
dens = mu/pi
print(f"(4) density mu/pi: std(log2) = {np.log2(dens).std():.4f}   "
      f"q01 = {np.quantile(dens, 0.01):.3f}   q99 = {np.quantile(dens, 0.99):.3f}")

phi = p_feed/rowsum

def theta_cascade(start, gens=10):
    w = start.copy()
    out = []
    for g in range(gens):
        th = float((w*phi).sum()/w.sum())
        out.append(th)
        wf = w*phi                     # mass taking the feed edge
        w2 = np.zeros(N)
        np.add.at(w2, sigma[has], wf[has])
        w = w2
        if w.sum() <= 0:
            break
    return out

th_mu = theta_cascade(mu)
th_pi = theta_cascade(pi)
print("(2) theta_g from mu:", "  ".join(f"{t:.4f}" for t in th_mu))
print("    theta_g from pi:", "  ".join(f"{t:.4f}" for t in th_pi))
print(f"    phi-bar = 1 - lam^-2 = {1-lam**-2:.4f}")
