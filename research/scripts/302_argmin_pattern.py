"""Classify the cb-argmin pattern: how many digits of j determine
pi(j) = argmin_e v[j + e*Nl]? Top digits vs bottom digits as predictors."""
import numpy as np
from math import log2

CACHE = "E:/projects/collatz/research/cache"

for lam in [1.05, 1.70]:
    for k in [12, 13, 14]:
        v = np.load(f"{CACHE}/v_lam{lam:.2f}_k{k}.npy")
        N = v.size; Nl = N//3
        cols = np.stack([v[:Nl], v[Nl:2*Nl], v[2*Nl:]], axis=1)
        pi = cols.argmin(axis=1)          # pattern in {0,1,2}, length Nl
        # tie check
        srt = np.sort(cols, axis=1)
        ties = (srt[:,1] - srt[:,0] < 1e-15).mean()
        counts = np.bincount(pi, minlength=3)/Nl

        # predictability: P(pi | m bottom digits) and P(pi | m top digits)
        j = np.arange(Nl, dtype=np.int64)
        nd = k-2  # ternary digits of j in [0, 3^(k-2))
        def acc_given(keyvals):
            # majority-class accuracy of pi given key
            order = np.argsort(keyvals, kind='stable')
            ks_ = keyvals[order]; ps = pi[order]
            acc = 0
            start = 0
            bounds = np.nonzero(np.diff(ks_))[0]+1
            segs = np.split(ps, bounds)
            for seg in segs:
                acc += np.bincount(seg, minlength=3).max()
            return acc/Nl
        line = []
        for m in [1, 2, 3, 4]:
            bot = j % 3**m
            top = j // 3**(nd-m)
            line.append(f"m={m}: bot={acc_given(bot):.3f} top={acc_given(top):.3f}")
        base = counts.max()
        print(f"lam={lam:.2f} k={k}: pi-dist={counts.round(3)} ties={ties:.1e} base={base:.3f} | " + " | ".join(line), flush=True)

        # self-similarity across k: compare pi at k with pi at k-1 (project j mod 3^(k-3))
        # (only when previous cached)
    print()
