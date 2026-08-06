"""
300_structure_battery.py
========================
Battery of structural hypotheses on cached eigenvectors (Script build_cache).
Digits are EXACTLY independent uniform under counting measure, so the ANOVA
decomposition is orthogonal and exact:
  Var(F) = sum_p main_p + sum_{p<q} I(p,q) + higher-order terms.

Tests per (lam, k):
 B1. Full pairwise interaction matrix I(p,q): banding — does I decay in
     distance d=q-p? Fit geometric rate; compare with 1/3 (carry-propagation
     of x4 in base 3) and with the prefix plateau.
 B2. Gibbs/Markov residual: total - [all mains + all adjacent pairs] —
     how close is F to a nearest-neighbour digit field?
 B3. Shape similarity: cosine similarity of main-effect profiles f_p across p.
 B4. Which component carries the prefix-ratio plateau: decay ratios of
     adjacent-pair interactions I(p,p+1) vs main_p vs plateau.
 B5. Noise floor: same statistics on digit-shuffled F.
 V.  Validation: class-1 columns = t x permuted class-0 (known exact).
"""
import numpy as np
from math import log2

ALPHA = log2(3.0)
CACHE = "E:/projects/collatz/research/cache"

def digits_of(N, k):
    i = np.arange(N, dtype=np.int64)
    D = np.empty((k-1, N), dtype=np.int8)
    x = i.copy()
    for p in range(k-1):
        D[p] = x % 3
        x //= 3
    return D

def anova(F, D, k):
    N = F.size
    mains = np.zeros(k-1)
    profiles = []
    for p in range(k-1):
        mm = np.array([F[D[p]==d].mean() for d in range(3)])
        profiles.append(mm - mm.mean())
        mains[p] = ((mm - mm.mean())**2).mean()
    # pairwise interactions (all pairs)
    P = k-1
    I = np.zeros((P, P))
    for p in range(P):
        for q in range(p+1, P):
            cell = np.zeros((3,3))
            for d1 in range(3):
                sel = D[p]==d1
                for d2 in range(3):
                    cell[d1,d2] = F[sel & (D[q]==d2)].mean()
            cell = cell - cell.mean(axis=1, keepdims=True) - cell.mean(axis=0, keepdims=True) + cell.mean()
            I[p,q] = (cell**2).mean()
    return mains, I, profiles

def report(tag, F, D, k, rho=None, lam=None):
    total = F.var()
    mains, I, prof = anova(F, D, k)
    P = k-1
    main_sum = mains.sum()
    adj = sum(I[p,p+1] for p in range(P-1))
    allpairs = I.sum()
    resid_markov = total - main_sum - adj
    resid_pairs = total - main_sum - allpairs
    print(f"\n=== {tag} ===")
    print(f"total Var = {total:.6f}")
    print(f"mains: {main_sum:.6f} ({100*main_sum/total:.2f}%)  "
          f"adj-pairs: {adj:.6f} ({100*adj/total:.2f}%)  "
          f"nonadj-pairs: {allpairs-adj:.6f} ({100*(allpairs-adj)/total:.2f}%)")
    print(f"Gibbs residual (beyond mains+adj): {resid_markov:.6f} ({100*resid_markov/total:.2f}%)")
    print(f"beyond-all-pairs (3-way+): {resid_pairs:.6f} ({100*resid_pairs/total:.2f}%)")
    # banding: mean I by distance
    print("I by distance d (mean over p):")
    prev = None
    for d in range(1, min(P, 7)):
        vals = [I[p,p+d] for p in range(P-d)]
        m = np.mean(vals)
        r = f" ratio={m/prev:.3f}" if prev and prev > 0 else ""
        print(f"  d={d}: {m:.3e}{r}")
        prev = m
    # adjacent-interaction decay in p
    adjseq = np.array([I[p,p+1] for p in range(P-1)])
    with np.errstate(divide='ignore', invalid='ignore'):
        adjrat = adjseq[1:]/adjseq[:-1]
    print("I(p,p+1) sequence:", " ".join(f"{x:.2e}" for x in adjseq))
    print("I(p,p+1) ratios:  ", " ".join(f"{r:.3f}" for r in adjrat))
    print("main_p ratios:    ", " ".join(f"{mains[p+1]/mains[p]:.3f}" if mains[p]>0 else "nan" for p in range(P-1)))
    # shape similarity of main profiles
    sims = []
    for p in range(P-1):
        a, b = prof[p], prof[p+1]
        na, nb = np.linalg.norm(a), np.linalg.norm(b)
        sims.append(float(a@b/(na*nb)) if na>0 and nb>0 else np.nan)
    print("main-profile cos(p,p+1):", " ".join(f"{s:+.3f}" for s in sims))

for lam in [1.05, 1.70]:
    k = 13
    fn = f"{CACHE}/v_lam{lam:.2f}_k{k}.npy"
    v = np.load(fn)
    rho = float(open(f"{CACHE}/rho_lam{lam:.2f}_k{k}.txt").read())
    N = v.size
    F = np.log2(v); F -= F.mean()
    D = digits_of(N, k)
    report(f"lam={lam} k={k}", F, D, k, rho, lam)

    # validation: class-1 = t * permuted class-0 (exact)
    t = (lam**-2)/rho
    Nl = N//3; Nl3 = Nl//3
    v0 = v[0::3]; v1 = v[1::3]
    s = np.arange(Nl)
    err = np.abs(v1 - t*v0[(4*s+2) % Nl]).max()
    print(f"validation class-1 copy: max err = {err:.2e} (should be ~1e-16)")

    # noise floor: digit-shuffled F
    rng = np.random.default_rng(42)
    Fs = F.copy(); rng.shuffle(Fs)
    mains_s, I_s, _ = anova(Fs, D, k)
    print(f"noise floor: main_p ~ {mains_s.mean():.2e}, I(p,q) ~ {I_s[I_s>0].mean():.2e}")
print("\nDONE")
