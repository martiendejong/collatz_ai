"""R323-329: STRUCTURE OF THE RESIDUAL FIELD (eigvec / roulette).
(1) digit-energy spectrum of the residual at k=13: complementary to roulette
    (coarse-digit dominated)?
(2) residual CV at k=13 vs k=15 vs k=17 at the same block depth mod 3^7:
    does the unexplained structure HOMOGENIZE with k? (direct view of theta)
(3) spatial autocorrelation of the residual: correlation length."""
import sys
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")

def theory_stationary(j):
    M = 3 ** j
    states = [r for r in range(M) if r % 3 != 0]
    idx = {r: i for i, r in enumerate(states)}
    P = np.zeros((len(states), len(states)))
    inv2 = pow((M + 1) // 2, 1, M)
    for r in states:
        b = (3 * r + 1) % M
        p = 0.5; x = b
        for w in range(1, 100):
            x = (x * inv2) % M
            P[idx[r], idx[x]] += p
            p *= 0.5
    P /= P.sum(1, keepdims=True)
    v = np.ones(len(states)) / len(states)
    for _ in range(4000): v = v @ P
    return dict(zip(states, v))

j = 7; M = 3 ** j
th = theory_stationary(j)
coset = np.array(sorted(s for s in th if s % 3 == 2))
t = np.array([th[s] for s in coset]); t /= t.sum()

print(f"{'k':>3} {'corr(log)':>9} {'CV(residual)':>12}")
res_by_k = {}
for k, path in ((13, "certificates/cert_k13.npy"), (15, "certificates/cert_k15.npy"),
                (17, "certificates/cert_k17.npy")):
    C = np.load(path, mmap_mode="r")
    N = 3 ** (k - 1)
    # block means by m mod 3^7: m = 3i+2 -> m mod M determined by i mod 3^{j-1}
    ev = np.asarray(C[:N], dtype=np.float64)
    ii = np.arange(N, dtype=np.int64) % (M // 3)
    sums = np.bincount(ii, weights=ev, minlength=M // 3)
    cnts = np.bincount(ii, minlength=M // 3)
    bm = sums / cnts   # indexed by i0 = (m mod M - 2)/3
    e = bm[(coset - 2) // 3]; e /= e.sum()
    res = e / t
    res_by_k[k] = res
    print(f"{k:>3} {np.corrcoef(np.log(e), np.log(t))[0,1]:>9.4f} {res.std()/res.mean():>12.4f}")

print("\nresidual digit-energy spectrum (k=13):")
L = np.log(res_by_k[13])
tot = L.var(); prev = tot; last = None
print(f"  total Var[log residual] = {tot:.5f}")
Lw = L.copy()
for d in range(j - 2):
    Lw = Lw.reshape(-1, 3).mean(1)
    v = Lw.var(); en = prev - v; prev = v
    r = f" ratio {en/last:.3f}" if last else ""
    print(f"  digit {d}: energy {en:.5f} ({en/tot*100:.1f}%){r}")
    last = en
print(f"  remaining coarse Var: {prev:.5f} ({prev/tot*100:.1f}%)")

x13 = res_by_k[13] - res_by_k[13].mean()
print("\nresidual lag autocorrelation (k=13):",
      " ".join(f"lag{l}:{(x13[:-l]*x13[l:]).mean()/(x13*x13).mean():+.3f}" for l in (1, 3, 9, 27, 81, 243)))
print("\ncross-k residual correlation (is the residual a stable object?):")
for a, b in ((13, 15), (15, 17), (13, 17)):
    r = np.corrcoef(res_by_k[a], res_by_k[b])[0, 1]
    print(f"  corr(res_k{a}, res_k{b}) = {r:.4f}")
