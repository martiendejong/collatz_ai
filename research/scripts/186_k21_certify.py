"""
186_k21_certify.py
==================
Final stage for k=21 (run after 185 converges with growth >= 1):
floor the converged memmap vector to integers (S=1e10) and verify ALL
3,486,784,401 constraints of L_21^NT(lambda_0) in exact integer arithmetic
with strict rational lower-bound weights (denominator 1e18, floor-minus-one)
-- the same protocol as 55e (k=20) and verify_certificates.py.

Usage: python 186_k21_certify.py [lam_num]      (default 1890 -> gamma 0.91836)
Reads research/k21/<cur>.npy per state.json; writes cert_k21.npy (int64
memmap, 28GB) next to it. Runtime estimate: 2-4h per billion constraints
(runbook) => ~7-14h, chunked, constant memory.
"""
import sys, os, json, math, decimal, time
import numpy as np

decimal.getcontext().prec = 80
K = 21
N = 3 ** (K - 1)
M = 3 ** K
Mc = 3 ** (K - 1)
Nl = N // 3
CH = 3 ** 12
S = 10 ** 10
Q = 10 ** 18

HERE = os.path.dirname(os.path.abspath(__file__))
DIR = os.path.abspath(os.path.join(HERE, "..", "k21"))
st = json.load(open(os.path.join(DIR, "state.json")))
lam_num = int(sys.argv[1]) if len(sys.argv) > 1 else 1890
lam_den = 1000

src = np.lib.format.open_memmap(os.path.join(DIR, st["cur"]), mode="r")
norm = st["norm"]
print(f"flooring {st['cur']} (norm {norm:.6f}) at S=1e10 ...", flush=True)
cert_path = os.path.join(DIR, "cert_k21.npy")
Ci = np.lib.format.open_memmap(cert_path, mode="w+", dtype=np.int64, shape=(N,))
for lo in range(0, N, CH * 27):
    hi = min(lo + CH * 27, N)
    Ci[lo:hi] = np.floor(np.asarray(src[lo:hi], dtype=np.float64) / norm * S).astype(np.int64)
Ci.flush()
del src
print("int certificate written; exact verify at "
      f"{lam_num}/{lam_den} (gamma = {math.log2(lam_num/lam_den):.5f})", flush=True)

dl = decimal.Decimal(lam_num) / decimal.Decimal(lam_den)
ln_l = dl.ln()
W0 = (lam_den ** 2 * Q) // (lam_num ** 2)
W2 = int(((decimal.Decimal('-0.41503749927884390') * ln_l).exp() * Q)
         .to_integral_value(rounding=decimal.ROUND_FLOOR)) - 1
W8 = int(((decimal.Decimal('0.58496250072115610') * ln_l).exp() * Q)
         .to_integral_value(rounding=decimal.ROUND_FLOOR)) - 1

C = np.lib.format.open_memmap(cert_path, mode="r")
viol = 0
t0 = time.time()
for lo in range(0, N, CH):
    hi = min(lo + CH, N)
    idx = np.arange(lo, hi, dtype=np.int64)
    m = 3 * idx + 2
    i4 = (((4 * m) % M) - 2) // 3
    mod9 = m % 9
    is2 = mod9 == 2
    is8 = mod9 == 8
    Cc = np.asarray(C[lo:hi]).astype(object)
    rhs = W0 * np.asarray(C[i4]).astype(object)
    for mask, mul, Wb in ((is2, 4, W2), (is8, 2, W8)):
        if mask.any():
            mm = m[mask]
            t = (((mul * mm - (2 if mul == 4 else 1)) // 3) % Mc)
            j = np.stack([(t - 2) // 3, ((t + Mc) - 2) // 3,
                          ((t + 2 * Mc) - 2) // 3])
            cb = np.minimum(np.minimum(np.asarray(C[j[0]]),
                                       np.asarray(C[j[1]])),
                            np.asarray(C[j[2]])).astype(object)
            rhs[mask] = rhs[mask] + Wb * cb
    viol += int((Cc * Q > rhs).sum())
    if (lo // CH) % 200 == 0:
        el = time.time() - t0
        print(f"  {hi:,}/{N:,} viol {viol}  [{el/60:.0f}m]", flush=True)

gamma = math.log2(lam_num / lam_den)
status = "VERIFIED *** NEW RECORD ***" if viol == 0 else "FAILED"
print(f"k=21 lam={lam_num}/{lam_den}: violations={viol}: {status} "
      f"pi(x) >= x^{gamma:.4f}", flush=True)
