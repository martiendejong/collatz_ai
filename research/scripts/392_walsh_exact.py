# 392: EXACT maximum linear correlation of the parity bit p_j over ALL linear
# masks, via fast Walsh-Hadamard on the full residue space mod 2^(j+1).
# p_j depends only on n mod 2^(j+1) (Terras), so enumerating all residues and
# transforming gives the EXACT max |corr| as a dyadic rational — upgrading the
# measured halving law (Obs 573/588) to certificate grade per depth.
import numpy as np

def fwht(a):
    h = 1
    n = len(a)
    while h < n:
        a = a.reshape(-1, 2, h)
        x = a[:, 0, :].copy()
        a[:, 0, :] = x + a[:, 1, :]
        a[:, 1, :] = x - a[:, 1, :]
        a = a.reshape(n)
        h *= 2
    return a

print(f"{'j':>3} {'max|corr| (alle maskers, odd seeds)':>36} {'als 2^-e':>9} {'wet 2^-ceil((j-3)/2)':>20}")
for j in range(1, 21):
    Nj = 1 << (j+1)
    x = np.arange(Nj, dtype=np.int64)
    n = x.copy()
    for _ in range(j):
        odd = (n & 1).astype(bool)
        n = np.where(odd, (3*n + 1) >> 1, n >> 1)
    p = (n & 1).astype(np.int64)
    # odd seeds only: function of bits 1..j
    podd = p[1::2]
    s = 1 - 2*podd.astype(np.float64)
    W = fwht(s.copy())
    mx = np.abs(W).max()/len(podd)
    # exact dyadic? mx * len = integer
    e = -np.log2(mx) if mx > 0 else np.inf
    wet = 2.0**(-max(0, -(-(j-3)//2)))
    print(f"{j:>3} {mx:>36.6f} {e:>9.3f} {wet:>20.4f}", flush=True)
