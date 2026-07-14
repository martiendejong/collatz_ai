"""R409: k=20 rerun with CHECKPOINTING (lesson from the crashed first run) +
chunked post-analysis + tempering exponent alpha_20. 60 iters (converged by 50)."""
import sys, math
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
ALPHA = math.log2(3)
k = 20
N = 3 ** (k - 1); M = 3 ** k; Mc = 3 ** (k - 1); M3 = 3 ** (k - 2)
LAM = 1.885
W0 = np.float32(LAM ** -2.0); W2 = np.float32(LAM ** (ALPHA - 2.0)); W8 = np.float32(LAM ** (ALPHA - 1.0))
CH = 3 ** 14
ITERS = 60
CKPT = "certificates/k20_eig.npy"

c = np.ones(N, dtype=np.float32)
print(f"k={k}: N={N:,}, lambda={LAM}, {ITERS} iters, checkpointing to {CKPT}", flush=True)
for it in range(ITERS):
    f = np.empty(N, dtype=np.float32)
    for lo in range(0, N, CH):
        hi = min(lo + CH, N)
        j = np.arange(lo, hi, dtype=np.int64)
        m = 3 * j + 2
        i4 = (((4 * m) % M) - 2) // 3
        out = W0 * c[i4]
        mod9 = m % 9
        for mask_val, mul, sub, W in ((2, 4, 2, W2), (8, 2, 1, W8)):
            sel = mod9 == mask_val
            if sel.any():
                mm = m[sel]
                t = ((mul * mm - sub) // 3) % Mc
                b = (t - 2) // 3
                cb = np.minimum(np.minimum(c[b], c[b + M3]), c[b + 2 * M3])
                out[sel] += W * cb
        f[lo:hi] = out
    del c
    mx = f.max(); f /= mx; c = f
    if it % 10 == 0:
        print(f"  iter {it}: norm {mx:.6f}", flush=True)
        np.save(CKPT, c)
np.save(CKPT, c)
print("eigenvector saved", flush=True)

# chunked q + CV_top
qs = 0.0; tot = 0.0; cv_acc = 0.0; cv_n = 0
for lo in range(0, M3, CH):
    hi = min(lo + CH, M3)
    t0, t1, t2 = c[lo:hi].astype(np.float64), c[lo+M3:hi+M3].astype(np.float64), c[lo+2*M3:hi+2*M3].astype(np.float64)
    tri = np.stack([t0, t1, t2])
    qs += tri.min(0).sum(); tot += tri.sum()
    mn = tri.mean(0); cv = tri.std(0) / mn
    cv_acc += cv.sum(); cv_n += cv.size
q = 3 * qs / tot; cv_top = cv_acc / cv_n
print(f"q = {q:.5f}   CV_top = {cv_top:.5f}", flush=True)

# CV profile per level p (chunked)
prof = []
for p in range(1, k - 1):
    B = 3 ** p
    acc_c = 0.0; acc_n = 0
    for lo in range(0, N, 3 * CH):
        hi = min(lo + 3 * CH, N)
        idx = np.arange(lo, hi, dtype=np.int64)
        sel = idx[(idx // B) % 3 == 0]
        sel = sel[sel + 2 * B < N]
        if sel.size == 0: continue
        t = np.stack([c[sel], c[sel + B], c[sel + 2 * B]]).astype(np.float64)
        cv = t.std(0) / t.mean(0)
        acc_c += cv.sum(); acc_n += cv.size
    prof.append(acc_c / acc_n)
    print(f"  profile p={p}: CV={prof[-1]:.4f}", flush=True)
prof = np.array(prof)
A = np.stack([prof[1:-3], prof[3:-1]], axis=1)
y = prof[2:-2]
coef, *_ = np.linalg.lstsq(A, y, rcond=None)
a, cc = coef
disc = 1 - 4 * a * cc
th = (1 - math.sqrt(disc)) / (2 * cc) if disc > 0 else float('nan')
print(f"\nFINAL k=20: q={q:.5f} CV_top={cv_top:.5f}")
print(f"lattice fit: a={a:.4f} c={cc:.4f} a+c={a+cc:.4f} theta={th:.4f}")
print("PREDICTIONS: (a,c)~(0.472,0.521), theta~0.849", flush=True)

# tempering exponent alpha_20: block means at mod 3^7 vs roulette
def theory_stationary(j):
    Mj = 3 ** j
    states = [r for r in range(Mj) if r % 3 != 0]
    idx = {r: i for i, r in enumerate(states)}
    P = np.zeros((len(states), len(states)))
    inv2 = pow((Mj + 1) // 2, 1, Mj)
    for r in states:
        b = (3 * r + 1) % Mj; pw = 0.5; x = b
        for w in range(1, 100):
            x = (x * inv2) % Mj
            P[idx[r], idx[x]] += pw; pw *= 0.5
    P /= P.sum(1, keepdims=True)
    v = np.ones(len(states)) / len(states)
    for _ in range(4000): v = v @ P
    return dict(zip(states, v))

jj = 7; Mj = 3 ** jj
th_d = theory_stationary(jj)
coset = np.array(sorted(s for s in th_d if s % 3 == 2))
tvec = np.array([th_d[s] for s in coset]); tvec /= tvec.sum()
B = Mj // 3
sums = np.zeros(B); cnts = np.zeros(B)
for lo in range(0, N, CH):
    hi = min(lo + CH, N)
    ii = np.arange(lo, hi, dtype=np.int64) % B
    sums += np.bincount(ii, weights=c[lo:hi].astype(np.float64), minlength=B)
    cnts += np.bincount(ii, minlength=B)
bm = sums / cnts
e = bm[(coset - 2) // 3]; e /= e.sum()
lt = np.log(tvec); le = np.log(e)
Amat = np.column_stack([lt - lt.mean(), np.ones_like(lt)])
(al, b0), *_ = np.linalg.lstsq(Amat, le, rcond=None)
r2 = 1 - ((le - Amat @ np.array([al, b0]))**2).sum()/((le - le.mean())**2).sum()
res = e / tvec
print(f"\nTEMPERING: alpha_20 = {al:.4f} (R2 {r2:.4f}); CV_res = {res.std()/res.mean():.4f}")
print("PREDICTION: alpha_20 ~ 0.878 (between k19 0.868 and k21), CV_res ~ 0.125", flush=True)
