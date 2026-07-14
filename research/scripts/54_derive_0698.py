"""R445-456: DERIVING THE 0.698 CONSTANT.
Min-Loss Identity (Thm 12): 1 = lam^-2 + (q/3)(lam^{a-2} + lam^{a-1}).
Implicit differentiation at the edge (lam=2, q=1) gives EXACT dgamma/dq.
Law candidate: (1 - gamma_k) = (dgamma/dq) * (1 - q_k). Measure q_k from all
four certificates (chunked) and test. Then 0.698 = 3.477*(1-q)/CV_res decomposes."""
import sys, math
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
A = math.log2(3)

# exact derivative at the edge
lam, q = 2.0, 1.0
F_q = (lam**(A-2) + lam**(A-1)) / 3
F_l = -2*lam**-3 + (q/3)*((A-2)*lam**(A-3) + (A-1)*lam**(A-2))
dl_dq = -F_q / F_l
dg_dq = dl_dq / (lam * math.log(2))
print(f"EXACT at the edge: dlam/dq = {dl_dq:.4f}, dgamma/dq = {dg_dq:.4f}")

print(f"\n{'k':>3} {'gamma':>7} {'1-gamma':>8} {'q':>7} {'1-q':>7} {'3.477*(1-q)':>11} {'ratio':>6}")
for k, lamk, path in ((13, 1.818, "certificates/cert_k13.npy"), (15, 1.841, "certificates/cert_k15.npy"),
                      (17, 1.86, "certificates/cert_k17.npy"), (19, 1.875, "certificates/cert_k19.npy")):
    C = np.load(path, mmap_mode="r")
    N = 3 ** (k - 1); M3 = N // 3
    CH = 3 ** 13
    qs = tot = 0.0
    for lo in range(0, M3, CH):
        hi = min(lo + CH, M3)
        t = np.stack([C[lo:hi], C[lo+M3:hi+M3], C[lo+2*M3:hi+2*M3]]).astype(np.float64)
        qs += t.min(0).sum(); tot += t.sum()
    qk = 3 * qs / tot
    g = math.log2(lamk)
    pred = dg_dq * (1 - qk)
    print(f"{k:>3} {g:>7.4f} {1-g:>8.4f} {qk:>7.4f} {1-qk:>7.4f} {pred:>11.4f} {(1-g)/pred:>6.3f}")
print("\nif ratio ~ 1: the gamma-lambda law's analytic skeleton = Min-Loss Identity derivative;")
print("0.698 = 3.477 x c1 x (CV_top/CV_res) decomposed.")
