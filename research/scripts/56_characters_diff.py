"""R510-521: (A) CHARACTER DECOMPOSITION of g. 2 is a primitive root mod 3^j, so
(Z/3^j)* = <2> cyclic; the coset r==2 mod 3 = odd discrete logs. Map g to dlog
coordinates and FFT: dominant frequencies identify g as a sum of characters.
(B) DIFFERENCE HUNT at roulette level: stationary law + digit-energy spectrum
for the 3n+1, 3n-1, 5n+1 roulettes -- is the roulette layer blind to the
map differences (expected), locating the true 3n+1-specific content?"""
import sys
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")

def roulette(j, mult, off, div=2):
    M = 3 ** j
    states = [r for r in range(M) if r % 3 != 0]
    idx = {r: i for i, r in enumerate(states)}
    P = np.zeros((len(states), len(states)))
    invd = pow(div, -1, M)
    for r in states:
        b = (mult * r + off) % M
        p = 0.5; x = b
        for w in range(1, 90):
            x = (x * invd) % M
            P[idx[r], idx[x]] += p; p *= 0.5
    P /= P.sum(1, keepdims=True)
    v = np.ones(len(states)) / len(states)
    for _ in range(4000): v = v @ P
    return states, v

j = 7; M = 3 ** j
states, v = roulette(j, 3, 1)
th = dict(zip(states, v))
coset = np.array(sorted(s for s in th if s % 3 == 2))
t = np.array([th[s] for s in coset]); t /= t.sum()
C = np.load("certificates/cert_k13.npy", mmap_mode="r")
N = 3 ** 12; B = M // 3
ii = np.arange(N, dtype=np.int64) % B
bm = np.bincount(ii, weights=np.asarray(C[:N], dtype=np.float64), minlength=B) / np.bincount(ii, minlength=B)
e = bm[(coset - 2) // 3]; e /= e.sum()
g = e / t; g = g / g.mean() - 1.0

# (A) dlog coordinates: dlog base 2 mod 3^j
order = 2 * 3 ** (j - 1)
dlog = {}
x = 1
for a in range(order):
    dlog[x] = a
    x = (x * 2) % M
# coset r==2 mod 3 <-> odd dlog; sequence over a = 1,3,5,...
seq = np.empty(order // 2)
for r in coset:
    seq[(dlog[r] - 1) // 2] = g[np.searchsorted(coset, r)]
F = np.fft.rfft(seq)
pw = np.abs(F) ** 2 / (np.abs(F) ** 2).sum()
top = np.argsort(-pw[1:])[:8] + 1
print("(A) character spectrum of g (fraction of power per frequency, dlog space):")
for f in top:
    print(f"  freq {f:>4}: {pw[f]:.4f}")
print(f"  top-3 frequencies carry {pw[top[:3]].sum():.3f} of total AC power; "
      f"total AC power fraction {(pw[1:].sum()):.4f}")

# (B) roulette layer difference hunt
print("\n(B) roulette stationary mod 9 + digit-energy for three maps:")
for name, mult, off in (("3n+1", 3, 1), ("3n-1", 3, -1), ("5n+1", 5, 1)):
    st9, v9 = roulette(2, mult, off)
    d = dict(zip(st9, v9))
    print(f"  {name}: mod9 = " + " ".join(f"{s}:{d[s]:.4f}" for s in sorted(d)))
    stj, vj = roulette(5, mult, off)
    cs = np.array(sorted(s for s in stj if s % 3 == 2))
    dd = dict(zip(stj, vj))
    L = np.log(np.array([dd[s] for s in cs]))
    tot = L.var(); Lw = L.copy(); prev = tot; ens = []
    for dgt in range(3):
        Lw = Lw.reshape(-1, 3).mean(1)
        vv = Lw.var(); ens.append(prev - vv); prev = vv
    print(f"       digit energies: {['%.4f' % x for x in ens]}  (Var {tot:.3f})")
