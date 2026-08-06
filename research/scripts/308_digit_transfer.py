"""
308_digit_transfer.py
=====================
The analytical cascade recursion, first-order digit transfer.

At the fixed point v = Fv/rho with F = A*P + CB (P = composition with T4,
CB = the cb-injection field). ANOVA is linear in the field, so EXACTLY:
    f_p(v) = [A * f_p(Pv) + f_p(CBv)] / rho          (digit-p main profiles)

We measure both contributions and then express f_p(Pv) in terms of the
f_q(v) profiles: the carry-mixing of x4+2 couples digit p of the image to
digits p, p-1 of the source. The resulting transfer operator's amplitude
eigenvalue should equal sqrt(prefix plateau): 0.64 (lam=1.05) / 0.837 (lam=1.70).
"""
import numpy as np
from math import log2

ALPHA = log2(3.0)
CACHE = "E:/projects/collatz/research/cache"

for lam in [1.05, 1.70]:
    k = 13
    v = np.load(f"{CACHE}/v_lam{lam:.2f}_k{k}.npy")
    rho = float(open(f"{CACHE}/rho_lam{lam:.2f}_k{k}.txt").read())
    A = lam**-2; B1 = lam**(ALPHA-2); B3 = lam**(ALPHA-1)
    N = v.size; Nl = N//3
    i = np.arange(N, dtype=np.int64)
    T4 = (4*i+2) % N
    s_arr, r_arr = np.divmod(i, 3)
    R1 = (4*s_arr) % Nl; R3 = (2*s_arr+1) % Nl
    cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])

    Pv = v[T4]
    CBv = np.zeros(N)
    CBv[r_arr == 0] = B1 * cb[R1[r_arr == 0]]
    CBv[r_arr == 2] = B3 * cb[R3[r_arr == 2]]
    # fixed point check: v = (A*Pv + CBv)/rho
    err = np.abs(A*Pv + CBv - rho*v).max() / v.mean()
    print(f"\n=== lam={lam} k={k}: fixed-point residual {err:.2e} ===")

    P = k-1
    digits = np.empty((P, N), dtype=np.int8)
    x = i.copy()
    for p in range(P):
        digits[p] = x % 3
        x //= 3

    def profiles(field):
        out = []
        for p in range(P):
            m = np.array([field[digits[p] == d].mean() for d in range(3)])
            out.append(m - m.mean())
        return out

    fv = profiles(v)
    fP = profiles(Pv)
    fC = profiles(CBv)

    print("p | ampl f_p(v) | A*ampl f_p(Pv)/rho | ampl f_p(CBv)/rho | check: A*fP+fC = rho*fv")
    for p in range(min(P, 6)):
        nv = np.linalg.norm(fv[p]); nP = np.linalg.norm(fP[p]); nC = np.linalg.norm(fC[p])
        lhs = A*fP[p] + fC[p]
        rhs = rho*np.array(fv[p])
        chk = np.abs(lhs - rhs).max() / (np.abs(rhs).max() + 1e-300)
        print(f"{p} | {nv:.5f} | {A*nP/rho:.5f} | {nC/rho:.5f} | rel err {chk:.1e}")

    # carry-mixing: project f_p(Pv) onto f_q(v) for q in {p-1, p, p+1}
    print("carry-mixing of P: f_p(Pv) decomposed on f_q(v) (cos overlaps):")
    for p in range(1, min(P, 5)):
        row = []
        for q in [p-1, p, p+1]:
            if 0 <= q < P:
                a = fP[p]; b = fv[q]
                na, nb = np.linalg.norm(a), np.linalg.norm(b)
                c = a@b/(na*nb) if na > 0 and nb > 0 else 0.0
                row.append(f"q={q}: cos={c:+.3f} amp={na/nb:.3f}")
        print(f"  p={p}: " + " | ".join(row))

    # empirical amplitude-transfer per digit: ratio of successive v-profile amplitudes
    amps = [np.linalg.norm(fv[p]) for p in range(P)]
    print("v-profile amplitude ratios f_{p+1}/f_p:",
          " ".join(f"{amps[p+1]/amps[p]:.3f}" if amps[p] > 1e-12 else "-" for p in range(min(P-1, 8))))
    print(f"sqrt(prefix plateau): {'0.640' if lam==1.05 else '0.837'}")
print("\nDONE")
