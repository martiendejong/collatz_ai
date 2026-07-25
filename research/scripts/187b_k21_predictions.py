"""
187b_k21_predictions.py
=======================
Test the remaining pre-registered k=21 predictions (PREDICTIONS.md, frozen
2026-07-16) on the converged k=21 vector (research/k21/va.npy, sweep 40):

  #2  alpha_21 = 0.887 +- 0.003   (tempering exponent, block mod 3^7)
  #3  CV_res(21) = 0.116 +- 0.004
  #4  theta(21) = 0.850 +- 0.001  (lattice fit, 25b windowed convention)
  #5  (a,c)(21) ~ (0.465, 0.528)
  (+ q and CV_top for the record; #7 already HIT: 0.97442 vs rule 0.97448)

Adapted from 25b_k20_resume.py post-analysis (same conventions), memmapped.
"""
import sys, math, json
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")

k = 21
N = 3 ** (k - 1)
M3 = 3 ** (k - 2)
CH = 3 ** 14

st = json.load(open("../k21/state.json"))
c = np.lib.format.open_memmap(f"../k21/{st['cur']}", mode="r")
print(f"k={k}: vector {st['cur']} (sweep {st['sweep']}), N={N:,}", flush=True)

# q + CV_top
qs = tot = cv_acc = 0.0
cv_n = 0
for lo in range(0, M3, CH):
    hi = min(lo + CH, M3)
    t0 = np.asarray(c[lo:hi], dtype=np.float64)
    t1 = np.asarray(c[lo + M3:hi + M3], dtype=np.float64)
    t2 = np.asarray(c[lo + 2 * M3:hi + 2 * M3], dtype=np.float64)
    tri = np.stack([t0, t1, t2])
    qs += tri.min(0).sum()
    tot += tri.sum()
    mn = tri.mean(0)
    cv = tri.std(0) / mn
    cv_acc += cv.sum()
    cv_n += cv.size
q = 3 * qs / tot
print(f"q = {q:.5f}   CV_top = {cv_acc/cv_n:.5f}", flush=True)

# CV profile per level p + lattice fit
prof = []
for p in range(1, k - 1):
    B = 3 ** p
    acc_c = 0.0
    acc_n = 0
    for lo in range(0, N, 3 * CH):
        hi = min(lo + 3 * CH, N)
        idx = np.arange(lo, hi, dtype=np.int64)
        sel = idx[(idx // B) % 3 == 0]
        sel = sel[sel + 2 * B < N]
        if sel.size == 0:
            continue
        t = np.stack([np.asarray(c[sel]), np.asarray(c[sel + B]),
                      np.asarray(c[sel + 2 * B])]).astype(np.float64)
        cv = t.std(0) / t.mean(0)
        acc_c += cv.sum()
        acc_n += cv.size
    prof.append(acc_c / acc_n)
    print(f"  profile p={p}: CV={prof[-1]:.4f}", flush=True)
prof = np.array(prof)
A = np.stack([prof[1:-3], prof[3:-1]], axis=1)
y = prof[2:-2]
coef, *_ = np.linalg.lstsq(A, y, rcond=None)
a, cc = coef
disc = 1 - 4 * a * cc
th = (1 - math.sqrt(disc)) / (2 * cc) if disc > 0 else float("nan")
print(f"lattice fit: a={a:.4f} c={cc:.4f} theta={th:.4f}")
print("PREDICTIONS #4/#5: theta = 0.850 +- 0.001, (a,c) ~ (0.465, 0.528)", flush=True)

# tempering alpha_21 (block means mod 3^7 vs roulette stationary)
def theory_stationary(j):
    Mj = 3 ** j
    states = [r for r in range(Mj) if r % 3 != 0]
    idx = {r: i for i, r in enumerate(states)}
    P = np.zeros((len(states), len(states)))
    inv2 = pow((Mj + 1) // 2, 1, Mj)
    for r in states:
        b = (3 * r + 1) % Mj
        pw = 0.5
        x = b
        for w in range(1, 100):
            x = (x * inv2) % Mj
            P[idx[r], idx[x]] += pw
            pw *= 0.5
    P /= P.sum(1, keepdims=True)
    v = np.ones(len(states)) / len(states)
    for _ in range(4000):
        v = v @ P
    return dict(zip(states, v))

jj = 7
Mj = 3 ** jj
th_d = theory_stationary(jj)
coset = np.array(sorted(s for s in th_d if s % 3 == 2))
tvec = np.array([th_d[s] for s in coset])
tvec /= tvec.sum()
B = Mj // 3
sums = np.zeros(B)
cnts = np.zeros(B)
for lo in range(0, N, CH):
    hi = min(lo + CH, N)
    ii = np.arange(lo, hi, dtype=np.int64) % B
    sums += np.bincount(ii, weights=np.asarray(c[lo:hi], dtype=np.float64),
                        minlength=B)
    cnts += np.bincount(ii, minlength=B)
bm = sums / cnts
e = bm[(coset - 2) // 3]
e /= e.sum()
lt = np.log(tvec)
le = np.log(e)
Amat = np.column_stack([lt - lt.mean(), np.ones_like(lt)])
(al, b0), *_ = np.linalg.lstsq(Amat, le, rcond=None)
r2 = 1 - ((le - Amat @ np.array([al, b0])) ** 2).sum() / ((le - le.mean()) ** 2).sum()
res = e / tvec
print(f"\nTEMPERING: alpha_21 = {al:.4f} (R2 {r2:.4f}); CV_res = {res.std()/res.mean():.4f}")
print("PREDICTIONS #2/#3: alpha_21 = 0.887 +- 0.003, CV_res = 0.116 +- 0.004", flush=True)
