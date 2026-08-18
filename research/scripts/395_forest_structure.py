# 395: dissect the self-similar forest (Obs 598): are the coherent masks at
# depth j+1 shifted copies of those at depth j (frontier translation)? And
# how does coherent energy distribute over mask weight classes?
import numpy as np
import math

def fwht_int(a):
    h = 1; n = len(a)
    while h < n:
        a = a.reshape(-1, 2, h)
        x = a[:, 0, :].copy()
        a[:, 0, :] = x + a[:, 1, :]
        a[:, 1, :] = x - a[:, 1, :]
        a = a.reshape(n); h *= 2
    return a

def coherent(j):
    Nj = 1 << (j+1)
    n = np.arange(Nj, dtype=np.int64)
    for _ in range(j):
        odd = (n & 1).astype(bool)
        n = np.where(odd, (3*n + 1) >> 1, n >> 1)
    s = (1 - 2*(n & 1)[1::2]).astype(np.int64)
    W = fwht_int(s)
    L = 1 << j
    corr2 = (W.astype(np.float64)/L)**2
    corr2[0] = 0
    thr = 16.0/L
    idx = np.nonzero(corr2 > thr)[0]
    return idx, corr2[idx]

sets = {}
for j in [19, 20, 21, 22]:
    idx, en = coherent(j)
    sets[j] = (set(int(x) for x in idx), {int(i): float(e) for i, e in zip(idx, en)})
    # gewichtsprofiel
    w = np.array([bin(int(i)).count('1') for i in idx])
    tot = en.sum()
    prof = {}
    for ww in range(int(w.min()), int(w.max())+1):
        m = w == ww
        if m.any():
            prof[ww] = float(en[m].sum())/tot
    top = sorted(prof.items(), key=lambda kv: -kv[1])[:5]
    print(f"j={j}: {len(idx)} coherent; energie-per-gewicht top: " +
          ", ".join(f"w={a}:{b:.2f}" for a, b in top), flush=True)

print()
for j in [19, 20, 21]:
    A, _ = sets[j]
    B, _ = sets[j+1]
    fr_same = len(B & A)/len(B)
    fr_shift = len(B & set(m << 1 for m in A))/len(B)
    fr_shift_or1 = len(B & (set(m << 1 for m in A) | set((m << 1) | 1 for m in A)))/len(B)
    both = len(B & (A | set(m << 1 for m in A) | set((m << 1) | 1 for m in A)))/len(B)
    print(f"j={j}->{j+1}: fractie van coherent(j+1) die = coherent(j): {fr_same:.3f} | "
          f"= shift(coherent(j)): {fr_shift:.3f} | shift+lowbit: {fr_shift_or1:.3f} | "
          f"verenigd: {both:.3f}", flush=True)
