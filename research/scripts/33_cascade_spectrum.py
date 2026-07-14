"""R259-265: THE CASCADE SPECTRUM. Two instruments on cert_k13/15/17:
(1) full per-digit correlation profile r(j) of the CV field (which 3-adic digit
    flip decorrelates how much) -- the ultrametric portrait;
(2) 'wavelet' energy spectrum of log-eigenvector: variance explained per digit
    position (coarse->fine), i.e. Var[ E[log C | digits >= j] ] increments.
    If the cascade is homogenizing (theta<1), fine-scale energy must decay
    geometrically; the decay ratio is an independent estimate of theta."""
import sys
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")

def analyze(k, path):
    C = np.load(path, mmap_mode="r")
    n = 3 ** (k - 1)
    # log field (float32 to save RAM at k=17)
    L = np.log(np.asarray(C[:n], dtype=np.float64))
    print(f"\n=== k={k}: {n:,} classes ===")
    # (2) energy spectrum: block-mean variances at successive 3-adic scales
    # index i corresponds to class m=3i+2; digit j of i = 3^j place.
    # coarse level j: average L over blocks where low j digits vary.
    tot_var = L.var()
    prev = 0.0
    print(f"  total Var[log C] = {tot_var:.6f}")
    print(f"  {'level j':>7} {'Var of block-means':>18} {'energy in digit j':>17} {'ratio':>7}")
    last_e = None
    energies = []
    Lw = L.copy()
    for j in range(0, k - 1):
        # block-average over lowest digit: reshape (n/3^{j+1}, 3) mean
        Lw = Lw.reshape(-1, 3).mean(1) if j > 0 else L.reshape(-1, 3).mean(1)
        v = Lw.var()
        e = (prev if j == 0 else prev) - v if j > 0 else tot_var - v
        # energy of digit j = Var at level j-1 minus Var at level j
        if j == 0:
            e = tot_var - v
        else:
            e = prev - v
        prev = v
        ratio = f"{e/last_e:.4f}" if last_e and last_e > 0 else "     -"
        energies.append(e)
        if j < 8 or j >= k - 4:
            print(f"  {j:>7} {v:>18.6f} {e:>17.6f} {ratio:>7}")
        last_e = e
    # geometric decay fit on fine half
    fine = [e for e in energies[:max(3,(k-1)//2)] if e > 0]
    if len(fine) >= 3:
        r = (fine[-1]/fine[0]) ** (1/(len(fine)-1))
        print(f"  fine-scale energy decay ratio (per digit) = {r:.4f}")
    return energies

e13 = analyze(13, "certificates/cert_k13.npy")
e15 = analyze(15, "certificates/cert_k15.npy")
import os
if os.path.exists("certificates/cert_k17.npy"):
    e17 = analyze(17, "certificates/cert_k17.npy")
