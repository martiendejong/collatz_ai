# 382: the discrete-log walk — the self-referential sampling in clock
# coordinates. Odd units mod 2^m = <-1> x <3> (ord(3) = 2^(m-2)); every odd n
# has clock position (eps, theta): n = (-1)^eps 3^theta mod 2^m. Pure x3 is
# theta -> theta+1 (rotation). Syracuse sampling makes theta JUMP.
# Measure: (1) equidistribution of theta along orbits; (2) independence of
# jump delta = theta' - theta - 1 from v; (3) jump distribution vs Haar;
# (4) sign-flip statistics. Any deviation = structure; none = camouflage
# extends to the clock coordinate.
import numpy as np
import random
from math import log2

M = 20
MOD = 1 << M
ORD = 1 << (M-2)
# dlog table for <3> mod 2^M
tab = {}
x = 1
for j in range(ORD):
    tab[x] = j
    x = (x*3) % MOD
def clock(n):
    r = n % MOD
    if r in tab:
        return 0, tab[r]
    return 1, tab[(-r) % MOD]

random.seed(382)
NORB, T = 400, 1600
thetas = []
deltas = []
vs = []
eps_flips = 0
tot = 0
for _ in range(NORB):
    n = random.getrandbits(900) | 1
    e0, th0 = clock(n)
    for t in range(T):
        m_ = 3*n + 1
        v = (m_ & -m_).bit_length() - 1
        n = m_ >> v
        e1, th1 = clock(n)
        thetas.append(th1)
        deltas.append((th1 - th0 - 1) % ORD)
        vs.append(v)
        eps_flips += (e1 != e0)
        e0, th0 = e1, th1
        tot += 1
thetas = np.array(thetas); deltas = np.array(deltas); vs = np.array(vs)

print(f"{tot} stappen, klok mod 2^{M} (ord = 2^{M-2})")
# (1) equidistribution of theta mod 2^j
print("\n(1) equidistributie theta mod 2^j (max afwijking van uniform, in sd-eenheden):")
for j in [2, 4, 6, 8]:
    cnt = np.bincount(thetas % (1 << j), minlength=1 << j)
    exp = tot/(1 << j)
    z = np.abs(cnt - exp).max()/np.sqrt(exp)
    print(f"   j={j}: max|dev| = {z:.2f} sd  ({'OK' if z < 4 else 'AFWIJKING'})")
# (2) jump vs v: mean jump per v-class mod small moduli
print("\n(2) sprong delta mod 8 per v-klasse (fracties; onafhankelijk => rijen gelijk):")
for v in [1, 2, 3, 4]:
    sel = deltas[vs == v] % 8
    fr = np.bincount(sel, minlength=8)/len(sel)
    print(f"   v={v}: {np.array2string(fr, precision=3)}  (n={len(sel)})")
# (3) jump distribution vs Haar at fine scale
print("\n(3) sprongverdeling mod 2^8 vs uniform:")
cnt = np.bincount(deltas % 256, minlength=256)
z = np.abs(cnt - tot/256).max()/np.sqrt(tot/256)
print(f"   max|dev| = {z:.2f} sd  ({'uniform' if z < 4 else 'STRUCTUUR'})")
# chi2-ish overall
chi = float(((cnt - tot/256)**2/(tot/256)).sum())
print(f"   chi2 = {chi:.0f} bij 255 vrijheidsgraden (verwacht ~255 +- 23)")
# (4) sign flips
print(f"\n(4) tekenwissel-fractie: {eps_flips/tot:.4f} (Haar: 0.5)")
# (5) autocorrelation of jumps (conspiracy would need persistent drift)
d_c = (deltas - deltas.mean())
ac1 = float((d_c[:-1]*d_c[1:]).mean()/(d_c**2).mean())
print(f"(5) autocorrelatie sprongen lag-1: {ac1:+.4f} (0 = geheugenloos)")
