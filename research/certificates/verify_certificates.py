"""Standalone verifier for the density-bound certificates (de Jong & Jengo, July 2026).

Verifies in EXACT INTEGER arithmetic that the deposited vectors are feasible solutions of the
Krasikov-Lagarias linear program L_k^NT(lambda_0) [KL2003, arXiv math/0205002, Sec. 2] at
  k=13, lambda_0 = 1818/1000  ==> pi(x) >= x^0.8624
  k=15, lambda_0 = 1841/1000  ==> pi(x) >= x^0.8805
via [KL2003, Theorem 2.2]. Requires only numpy and Python >= 3.8. Runtime: ~2 minutes.

Weights are replaced by strict rational LOWER bounds (denominator 10^18) computed with 80-digit
decimal arithmetic and floor-minus-one; alpha = log2(3) bracketed in
(1.5849625007211561, 1.5849625007211563). Any feasible c with the smaller weights is feasible
with the true weights, so the verification errs strictly against the claim.
"""
import math, decimal
import numpy as np

decimal.getcontext().prec = 80

def verify(k, lam_num, lam_den, path, chunk=3 ** 12):
    """Chunked exact-integer verification (constant memory beyond the certificate itself)."""
    N = 3 ** (k - 1); M = 3 ** k; Mc = 3 ** (k - 1)
    Q = 10 ** 18
    dl = decimal.Decimal(lam_num) / decimal.Decimal(lam_den)
    ln_l = dl.ln()
    W0 = (lam_den ** 2 * Q) // (lam_num ** 2)
    e2 = decimal.Decimal('-0.41503749927884390')   # < alpha - 2  (more negative)
    e8 = decimal.Decimal('0.58496250072115610')    # < alpha - 1
    W2 = int(((e2 * ln_l).exp() * Q).to_integral_value(rounding=decimal.ROUND_FLOOR)) - 1
    W8 = int(((e8 * ln_l).exp() * Q).to_integral_value(rounding=decimal.ROUND_FLOOR)) - 1
    C = np.load(path, mmap_mode='r')
    viol = 0
    for lo in range(0, N, chunk):
        hi = min(lo + chunk, N)
        idx = np.arange(lo, hi, dtype=np.int64)
        m = 3 * idx + 2
        i4 = (((4 * m) % M) - 2) // 3
        mod9 = m % 9
        is2 = mod9 == 2; is8 = mod9 == 8
        Cc = np.asarray(C[lo:hi]).astype(object)
        rhs = W0 * np.asarray(C[i4]).astype(object)
        for mask, mul, Wb in ((is2, 4, W2), (is8, 2, W8)):
            if mask.any():
                mm = m[mask]
                t = (((mul * mm - (2 if mul == 4 else 1)) // 3) % Mc)
                j = np.stack([(t - 2) // 3, ((t + Mc) - 2) // 3, ((t + 2 * Mc) - 2) // 3])
                cb = np.minimum(np.minimum(np.asarray(C[j[0]]), np.asarray(C[j[1]])),
                                np.asarray(C[j[2]])).astype(object)
                rhs[mask] = rhs[mask] + Wb * cb
        viol += int((Cc * Q > rhs).sum())
    gamma = math.log2(lam_num / lam_den)
    status = "VERIFIED" if viol == 0 else "FAILED"
    print(f"k={k}, lambda={lam_num}/{lam_den}: {N:,} constraints, violations={viol}: {status}"
          f" => pi(x) >= x^{gamma:.4f}")
    return viol == 0

if __name__ == "__main__":
    import os
    results = [verify(13, 1818, 1000, "cert_k13.npy"),
               verify(15, 1841, 1000, "cert_k15.npy")]
    if os.path.exists("cert_k17.npy"):
        results.append(verify(17, 186, 100, "cert_k17.npy"))
    if os.path.exists("cert_k19.npy"):
        results.append(verify(19, 15, 8, "cert_k19.npy"))
    print("ALL CERTIFICATES VALID" if all(results) else "VERIFICATION FAILURE")
