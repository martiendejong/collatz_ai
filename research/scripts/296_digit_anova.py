import numpy as np
from math import log2

ALPHA = log2(3.0)

for lam in [1.05, 1.70]:
    A = lam**-2; B1 = lam**(ALPHA-2); B3 = lam**(ALPHA-1)
    k = 12
    N = 3**(k-1); Nl = N//3
    i = np.arange(N, dtype=np.int64)
    T4 = (4*i+2) % N
    s_arr, r_arr = np.divmod(i, 3)
    m0 = r_arr==0; m2 = r_arr==2
    R1 = (4*s_arr) % Nl; R3 = (2*s_arr+1) % Nl
    del s_arr, r_arr
    v = np.ones(N)
    for _ in range(500):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        w_ = A*v[T4]
        w_[m2] += B3*cb[R3[m2]]
        w_[m0] += B1*cb[R1[m0]]
        v = w_ / w_.max()
    F = np.log2(v)
    F = F - F.mean()
    total = F.var()

    # ternary digits of index i: digit 0 = least significant
    digits = np.empty((k-1, N), dtype=np.int8)
    x = i.copy()
    for p in range(k-1):
        digits[p] = x % 3
        x //= 3

    # main effects
    main = np.zeros(k-1)
    effects = []
    resid = F.copy()
    for p in range(k-1):
        means = np.array([F[digits[p]==d].mean() for d in range(3)])
        effects.append(means)
        main[p] = (means**2).mean()   # variance of the main effect (equal class sizes)
    main_sum = main.sum()

    # pairwise interactions for adjacent digit pairs (p, p+1)
    inter = np.zeros(k-2)
    for p in range(k-2):
        cell = np.zeros((3,3))
        for d1 in range(3):
            for d2 in range(3):
                cell[d1,d2] = F[(digits[p]==d1)&(digits[p+1]==d2)].mean()
        cell = cell - cell.mean(axis=1, keepdims=True) - cell.mean(axis=0, keepdims=True) + cell.mean()
        inter[p] = (cell**2).mean()
    inter_sum = inter.sum()

    print(f"lam={lam}: total Var(log2 v) = {total:.5f}")
    print(f"  main effects sum      = {main_sum:.5f} ({100*main_sum/total:.1f}%)")
    print(f"  adjacent-pair inter.  = {inter_sum:.5f} ({100*inter_sum/total:.1f}%)")
    print(f"  per-digit main effect (p=0 least significant):")
    print("   ", " ".join(f"{m:.4f}" for m in main))
    print(f"  main-effect decay ratios (p->p+1):", 
          " ".join(f"{main[p+1]/main[p]:.3f}" for p in range(k-2)))
    # digit-profile shapes: are the 3-level means similar across positions (self-similar)?
    for p in [2, 5, 8]:
        e = effects[p]
        print(f"  digit {p} means: {e[0]:+.4f} {e[1]:+.4f} {e[2]:+.4f}")
