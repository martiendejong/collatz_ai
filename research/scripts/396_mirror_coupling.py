# 396: the MIRROR COUPLING A_j = <s_j(n) s_j(-n)> — mechanism for odd-weight
# dominance. Complementing bits 1..j of an odd residue is negation mod 2^(j+1),
# and for Boolean spectra: E_odd - E_even = -A_j. Prediction: A_1 = -1 exactly
# (the +/- pair: exactly one of (3n+1)/2, (3n-1)/2 is odd), then decay.
import numpy as np

print(f"{'j':>3} {'A_j (spiegel)':>13} {'oneven-fractie':>14} {'|A|-ratio':>9}")
prev = None
for j in range(1, 25):
    Nj = 1 << (j+1)
    n = np.arange(Nj, dtype=np.int64)
    for _ in range(j):
        odd = (n & 1).astype(bool)
        n = np.where(odd, (3*n + 1) >> 1, n >> 1)
    p = (n & 1).astype(np.int8)
    s = 1 - 2*p.astype(np.float64)
    idx = np.arange(Nj)
    neg = (Nj - idx) % Nj
    m = (idx & 1).astype(bool)      # oneven residuen
    A = float((s[m]*s[neg[m]]).mean())
    frac_odd = (1 - A)/2
    r = f"{abs(A)/abs(prev):.3f}" if prev not in (None, 0.0) else "  -"
    print(f"{j:>3} {A:>13.6f} {frac_odd:>14.6f} {r:>9}", flush=True)
    prev = A
