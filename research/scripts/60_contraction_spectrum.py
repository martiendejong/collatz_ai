"""R551-560: THE CONTRACTION SPECTRUM. Build the LINEARIZED K-L operator at the
edge (lam=2: W0=1/4, W2=3/4, W8=3/2; min -> triple mean) on classes mod 3^j,
compute its full spectrum. Leading eigenvalue ~1 = the roulette direction;
the SUBDOMINANT eigenvalue = the linear contraction rate of the g-modes.
Compare with the measured lambda-ratio (~0.93/digit, ~0.87/2digits)."""
import sys, math
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
ALPHA = math.log2(3)
lam = 2.0
W0 = lam**-2; W2 = lam**(ALPHA-2); W8 = lam**(ALPHA-1)

for j in (4, 5, 6):
    M = 3 ** j; Mc = 3 ** (j - 1)
    cls = [m for m in range(M) if m % 3 == 2]
    idx = {m: i for i, m in enumerate(cls)}
    n = len(cls)
    L = np.zeros((n, n))
    for m in cls:
        i = idx[m]
        # W0 term: class 4m mod 3^j
        L[i, idx[(4 * m) % M if (4 * m) % M % 3 == 2 else ((4 * m) % M)]] += W0
        # branch term for m == 2 or 8 mod 9: triple mean over children
        r9 = m % 9
        if r9 == 2:
            t = ((4 * m - 2) // 3) % Mc
            W = W2
        elif r9 == 8:
            t = ((2 * m - 1) // 3) % Mc
            W = W8
        else:
            continue
        for a in range(3):
            ch = (t + a * Mc) % M
            L[i, idx[ch]] += W / 3
    ev = np.linalg.eigvals(L)
    ev = ev[np.argsort(-np.abs(ev))]
    print(f"j={j} ({n} classes): |eigenvalues| top 6: " +
          " ".join(f"{abs(e):.4f}" for e in ev[:6]))
    print(f"   leading {abs(ev[0]):.4f} (roulette direction); "
          f"SUBDOMINANT = {abs(ev[1]):.4f}  <- linear contraction rate")
print("\nmeasured lambda-ratio: ~0.93/digit (0.87/2digits); creep 0.87->0.89")
