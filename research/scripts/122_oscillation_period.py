"""
122_oscillation_period.py
==========================
Measure the oscillation period of the dominant complex eigenmode
of the Collatz macro-step chain at multiple moduli.

Key question: is the oscillation period (~4-5 macro-steps from mod-256 analysis)
universal across moduli, or does it scale with the modulus?
"""
import numpy as np, math
from scipy.sparse.linalg import eigs
import scipy.sparse as sp

def v2(x):
    if x == 0: return 999
    c = 0
    while x % 2 == 0: x //= 2; c += 1
    return c

def macro_step(n):
    K = v2(n+1); m = (n+1) >> K; x = m*(3**K)-1; l = v2(x)
    return x >> l, K, l

print("Oscillation period of dominant complex eigenmode:")
print(f"{'Mod':>8} {'N':>6} {'Re(l2)':>10} {'Im(l2)':>10} {'|l2|':>8} {'phase_deg':>12} {'period':>10}")
print("-" * 75)

configs = [(256, 8192), (512, 4096), (1024, 2048), (2048, 1024)]
for MOD, NS in configs:
    odd_res = list(range(1, MOD, 2)); N = len(odd_res)
    idx = {r: i for i, r in enumerate(odd_res)}
    P = np.zeros((N, N), dtype=np.float32)
    for i, r in enumerate(odd_res):
        K0 = v2(r+1); counts = np.zeros(N, dtype=np.float32); valid = 0
        for k in range(NS):
            n = r + MOD*k
            if v2(n+1) != K0: continue
            n_out, _, _ = macro_step(n)
            r_out = n_out % MOD
            if r_out % 2 == 1 and r_out in idx:
                counts[idx[r_out]] += 1; valid += 1
        if valid > 0: P[i] = counts / valid
    P64 = P.astype(np.float64)
    P_sp = sp.csr_matrix(P64.T)
    vals, _ = eigs(P_sp, k=10, which="LM")
    vals_by_abs = sorted(vals, key=lambda x: -abs(x))
    ev = vals_by_abs[1]
    rho = abs(ev)
    theta = abs(math.atan2(ev.imag, ev.real))
    period = 2 * math.pi / theta if theta > 0.001 else float('inf')
    print(f"{MOD:>8} {N:>6} {ev.real:>10.5f} {ev.imag:>10.5f} {rho:>8.5f} {math.degrees(theta):>12.2f} {period:>10.3f}")
    # Also show top 5 complex modes
    print(f"         Top 5 eigenvalues by |lambda|:")
    for ev2 in vals_by_abs[:6]:
        rho2 = abs(ev2)
        if rho2 < 0.001: continue
        theta2 = abs(math.atan2(ev2.imag, ev2.real))
        period2 = 2*math.pi/theta2 if theta2 > 0.001 else float('inf')
        print(f"           lambda = {ev2.real:+.5f} + {ev2.imag:+.5f}i  "
              f"|l|={rho2:.5f}  period={period2:.2f} steps")
    print()

# Also: K-autocorrelation to see if 4-5 steps is the correlation length
print("=" * 70)
print("K-VALUE AUTOCORRELATION (lag 1..10)")
print("=" * 70)
MOD = 256
NS_LONG = 1000  # long orbit segments
n_start = 7  # start from an odd number
orbit_K = []
n = n_start
for _ in range(NS_LONG * 200):
    K, _, _ = v2(n+1), None, None
    K = v2(n+1)
    orbit_K.append(K)
    n_out, _, _ = macro_step(n)
    n = n_out
    if n == 1: break

# Autocorrelation
K_arr = np.array(orbit_K[:1000], dtype=float)
K_mean = K_arr.mean()
K_var = K_arr.var()
print(f"Mean K = {K_mean:.4f} (theoretical: 2.0)")
print(f"Var K = {K_var:.4f}")
print()
print(f"{'Lag':>6} {'Autocorr':>12}")
for lag in range(1, 12):
    if lag >= len(K_arr): break
    cov = np.cov(K_arr[:-lag], K_arr[lag:])[0, 1]
    acf = cov / K_var
    print(f"{lag:>6} {acf:>12.5f}")
