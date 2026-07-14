"""R302-305: THE EIGEN-VISIT BRIDGE at mod 27.
R284-289 found the K-L eigenvector's mod-9 block means track the forward
orbit-visit stationary law within 10%. Sharpen: measure both fields on all
9 classes mod 27 (m == 2 mod 3) and correlate. If the bridge is real, the
correlation should persist at the finer digit."""
import sys, random
import numpy as np
from collections import Counter
sys.stdout.reconfigure(encoding="utf-8")
random.seed(302)

# forward stationary law mod 27
stat = Counter()
for _ in range(80_000):
    n = random.randrange(3, 1 << 44) | 1
    m = n
    while m != 1:
        m = 3*m+1 if m % 2 else m//2
        if m % 2 and m != 1: stat[m % 27] += 1
tot = sum(stat.values())

C = np.load("certificates/cert_k13.npy", mmap_mode="r")
k = 13; n = 3 ** (k - 1)
idx = np.arange(n, dtype=np.int64)
m = 3 * idx + 2
mod27 = m % 27
vals = np.asarray(C[:n], dtype=np.float64)

classes = [c for c in range(27) if c % 3 == 2]
ev, fw = [], []
print(f"{'mod27':>5} {'eigvec share':>12} {'fwd-visit share':>15}")
evm = {c: vals[mod27 == c].mean() for c in classes}
evtot = sum(evm.values())
fwtot = sum(stat[c] for c in classes)
for c in classes:
    e = evm[c]/evtot; f = stat[c]/fwtot
    ev.append(e); fw.append(f)
    print(f"{c:>5} {e:>12.4f} {f:>15.4f}")
r = np.corrcoef(ev, fw)[0, 1]
rl = np.corrcoef(np.log(ev), np.log(fw))[0, 1]
print(f"\ncorrelation: linear r = {r:.4f}, log-log r = {rl:.4f}  (9 classes mod 27)")
