# 393: which masks WIN the Walsh maximum? Structure of the argmax masks of
# the parity function per depth j, plus deeper exact maxima (int FWHT).
import numpy as np

def fwht_int(a):
    h = 1; n = len(a)
    while h < n:
        a = a.reshape(-1, 2, h)
        x = a[:, 0, :].copy()
        a[:, 0, :] = x + a[:, 1, :]
        a[:, 1, :] = x - a[:, 1, :]
        a = a.reshape(n); h *= 2
    return a

print(f"{'j':>3} {'max|corr|':>10} {'winnend masker (bits 1..j)':>28} {'gewicht':>7} {'runs':>5}")
rates = {}
for j in range(6, 27):
    Nj = 1 << (j+1)
    n = np.arange(Nj, dtype=np.int64)
    for _ in range(j):
        odd = (n & 1).astype(bool)
        n = np.where(odd, (3*n + 1) >> 1, n >> 1)
    s = (1 - 2*(n & 1)[1::2]).astype(np.int64 if j < 24 else np.int32)
    del n
    W = fwht_int(s)
    del s
    L = 1 << j
    am = int(np.abs(W).argmax())
    mx = abs(int(W[am]))/L
    rates[j] = mx
    mb = bin(am)[2:].zfill(j)[::-1]  # bit i = seed bit i+1
    wt = mb.count('1')
    runs = len([r for r in mb.split('0') if r])
    show = mb if j <= 26 else mb[:26]
    print(f"{j:>3} {mx:>10.6f} {show:>28} {wt:>7} {runs:>5}", flush=True)
    del W
ks = sorted(rates)
import math
r = [rates[ks[i+1]]/rates[ks[i]] for i in range(len(ks)-1)]
tail = r[-8:]
gm = math.exp(sum(math.log(x) for x in tail)/len(tail))
print(f"\nstaart-vervalvoet (laatste 8 ratio's, geometrisch): {gm:.4f}")
