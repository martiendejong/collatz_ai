"""R595-600: THE 2-ADIC MIRROR. The survivor indicator field on odd residues
mod 2^K: digit-energy spectrum (2-adic wavelet decomposition) + conditional
survivor probability per bit. Does the 2-adic side carry a cascade like the
3-adic eigenvector, and where does survivor information live?"""
import sys, math
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")

def survivors(K):
    cur = [(1, 0, 0, 0)]
    for k in range(1, K + 1):
        out = []
        M = 1 << k
        for r, u, j, b in cur:
            for r2 in (r, r + (1 << (k - 1))):
                uu, jj, bb = u, j, b
                ok = True
                while jj < k and ok:
                    val = (pow(3, uu, M << 1) * r2 + bb) >> jj
                    if val & 1:
                        uu += 1; bb = 3 * bb + (1 << jj)
                    jj += 1
                    if 3 ** uu < (1 << jj): ok = False
                if ok: out.append((r2, uu, jj, bb))
        cur = out
    return set(r for r, *_ in cur)

K = 18
S = survivors(K)
n_odd = 1 << (K - 1)
print(f"survivors mod 2^{K}: {len(S):,} of {n_odd:,} odd (density {len(S)/n_odd:.5f})")

# indicator field over odd residues, indexed by i = (r-1)/2
F = np.zeros(n_odd)
for r in S: F[(r - 1) >> 1] = 1.0
tot = F.var()
print(f"digit-energy spectrum of survivor indicator (total Var {tot:.5f}):")
prev = tot; Fw = F.copy(); last = None
for d in range(K - 6):
    Fw = Fw.reshape(-1, 2).mean(1)
    v = Fw.var(); e = prev - v; prev = v
    r = f" ratio {e/last:.3f}" if last else ""
    print(f"  bit {d+1}: energy {e:.5f} ({e/tot*100:.1f}%){r}")
    last = e
print(f"  remaining coarse Var: {prev:.5f} ({prev/tot*100:.1f}%)")

# survivor probability by single bit value (which bits carry information?)
print("\nP(survive | bit b of r) informativeness:")
idx = np.arange(n_odd)
r = 2 * idx + 1
for b in range(1, 10):
    bit = (r >> b) & 1
    p1 = F[bit == 1].mean(); p0 = F[bit == 0].mean()
    print(f"  bit {b}: P(surv|1)={p1:.4f} P(surv|0)={p0:.4f}  ratio {p1/max(p0,1e-9):.2f}")
