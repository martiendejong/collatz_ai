"""R316-322: THE RESIDUAL FIELD. Decompose eigenvector = roulette x residual.
(1) compute roulette stationary at mod 3^7 (1458 states);
(2) eigenvector block means at mod 3^7 -> residual = eig/theory;
(3) CV of raw vs residual block means: how much structure does the roulette explain?
(4) digit-energy spectrum of the ROULETTE measure itself: does it reproduce
    the 0.20/digit cascade ratio of R259-265? If yes, the cascade is roulette,
    and only the residual carries the open spectral content."""
import sys
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")

def theory_stationary(j):
    M = 3 ** j
    states = [r for r in range(M) if r % 3 != 0]
    idx = {r: i for i, r in enumerate(states)}
    n = len(states)
    P = np.zeros((n, n))
    inv2 = pow((M + 1) // 2, 1, M)
    for r in states:
        b = (3 * r + 1) % M
        p = 0.5; x = b
        for w in range(1, 100):
            x = (x * inv2) % M
            P[idx[r], idx[x]] += p
            p *= 0.5
    P /= P.sum(1, keepdims=True)
    v = np.ones(n) / n
    for _ in range(4000): v = v @ P
    return states, v

j = 7; M = 3 ** j
states, v = theory_stationary(j)
th = dict(zip(states, v))

C = np.load("certificates/cert_k13.npy", mmap_mode="r")
k = 13; N = 3 ** (k - 1)
mvals = 3 * np.arange(N, dtype=np.int64) + 2
ev = np.asarray(C[:N], dtype=np.float64)

cls = [s for s in states if s % 3 == 2]
mm = mvals % M
e = np.array([ev[mm == s].mean() for s in cls])
t = np.array([th[s] for s in cls])
e /= e.sum(); t /= t.sum()
res = e / t
print(f"mod 3^{j} = {M}: {len(cls)} classes ==2 mod 3")
print(f"  corr(eig, theory) = {np.corrcoef(e, t)[0,1]:.4f} (log-log {np.corrcoef(np.log(e), np.log(t))[0,1]:.4f})")
print(f"  CV of eig block means      = {e.std()/e.mean():.4f}")
print(f"  CV of theory               = {t.std()/t.mean():.4f}")
print(f"  CV of RESIDUAL (eig/theory)= {res.std()/res.mean():.4f}  <- unexplained structure")

# (4) digit-energy spectrum of the roulette measure (log field over classes mod 3^7)
L = np.log(np.array([th[s] for s in sorted(states)]))
# order states by 3-adic index: sorted(states) is by value; digits of (s) low->high
# arrange into array indexed by (s - offset)/? -- states are all r % 3 != 0: two cosets.
# use coset r==2 mod 3 (matches cert indexing): index i = (r-2)/3, i in 0..3^{j-1}-1
coset = np.array(sorted(s for s in states if s % 3 == 2))
Lc = np.log(np.array([th[s] for s in coset]))
tot = Lc.var()
print(f"\nroulette log-measure digit-energy (coset ==2 mod 3), total Var = {tot:.4f}:")
Lw = Lc.copy(); prev = tot; last = None
for d in range(j - 1):
    Lw = Lw.reshape(-1, 3).mean(1)
    vvar = Lw.var()
    en = prev - vvar; prev = vvar
    ratio = f"  ratio {en/last:.4f}" if last else ""
    print(f"  digit {d}: energy {en:.5f}{ratio}")
    last = en
print("\n[compare: eigenvector cascade had ~89% in digit 0 and decay ratio ~0.20/digit]")
