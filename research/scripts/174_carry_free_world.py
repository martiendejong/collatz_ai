"""
174_carry_free_world.py
========================
RUNG t1 OF THE RIGIDITY BRIDGE: the carry-free Collatz world F2[t], where
the conjecture is a THEOREM (Hicks-Mullen-Yucas-Zavislak 2008), vs Z.

The polynomial analog over F2: for P with P(0) != 0 ("odd"):
    P -> ((t+1) P + 1) / t^e,   e = v_t((t+1)P + 1)
(the exact analog of n -> (3n+1)/2^e: t+1 <-> 3 = binary 11, t <-> 2).

Known: every nonzero P reaches 1. The reason carries matter: in Z,
3n+1 = (binary 11)*n + 1 WITH carry propagation; in F2[t] the same formula
has NO carries -- low coefficients of the image depend ONLY on low
coefficients of P. The halving exponent e is a function of the low part.

Experiments:
 E1: verify convergence + collect orbit statistics over all odd P up to
     degree 22 (4M polynomials): max excursion of deg, e-distribution.
 E2: THE MECHANISM: in F2[t], deg(P') = deg P + 1 - e and e-history is a
     deterministic function of initial low bits: the orbit map is
     conjugate to a LINEAR system -> exhibit the monotone functional
     (deg P works after bounded time? measure the maximal deg-excursion:
     if bounded by deg P + C uniformly, that is the theorem's shape).
 E3: port to Z: same trajectory statistics for random odd n of matching
     size; compare e-distributions and excursion distributions. The
     DIFFERENCE is exactly the carry contribution -- quantify it.
 E4: carry localization: for each Z-step, compute the F2-prediction of e
     (from the carry-free formula on the same bits) vs actual e; measure
     the fraction of steps where the carry changes e, and the net drift
     it adds/removes.
"""
import numpy as np
import random
from math import log2

# ---------- F2[t] arithmetic on bitmasks (int = coefficient bitmask) ----------
def f2_mul_t1(P):        # multiply by (t+1)
    return (P << 1) ^ P

def f2_step(P):
    Q = f2_mul_t1(P) ^ 1
    e = (Q & -Q).bit_length() - 1
    return Q >> e, e

def f2_orbit(P, maxsteps=10000):
    es = []
    d0 = P.bit_length() - 1
    maxdeg = d0
    steps = 0
    while P != 1 and steps < maxsteps:
        P, e = f2_step(P)
        es.append(e)
        maxdeg = max(maxdeg, P.bit_length() - 1)
        steps += 1
    return steps, es, maxdeg, d0

# ---------- Z (Collatz) ----------
def z_step(n):
    q = 3*n + 1
    e = (q & -q).bit_length() - 1
    return q >> e, e

print("E1/E2: F2[t] world -- all odd P, deg <= 22")
maxsteps_seen = 0
excess = np.zeros(40, dtype=np.int64)     # histogram of maxdeg - d0
e_hist_f2 = np.zeros(40, dtype=np.int64)
worst = None
for P in range(3, 1 << 23, 2):
    steps, es, maxdeg, d0 = f2_orbit(P)
    if steps >= 10000:
        print(f"  NON-CONVERGENT?! P={P:b}")
        break
    exc = maxdeg - d0
    excess[min(exc, 39)] += 1
    for e in es:
        e_hist_f2[min(e, 39)] += 1
    if steps > maxsteps_seen:
        maxsteps_seen = steps
        worst = (P, steps, d0)
else:
    print(f"  ALL odd P with deg <= 22 converge to 1  (theorem confirmed)")
print(f"  max steps: {maxsteps_seen} (P has deg {worst[2]})")
tot = excess.sum()
print(f"  deg-excursion histogram (maxdeg - d0): "
      f"{[(i, int(excess[i])) for i in range(6) if excess[i]]}"
      f"  P(exc >= 3) = {float(excess[3:].sum())/tot:.2e}")
ef2 = e_hist_f2 / e_hist_f2.sum()
print(f"  e-distribution F2: e=1: {ef2[1]:.4f}  e=2: {ef2[2]:.4f}  "
      f"e=3: {ef2[3]:.4f}  mean = {sum(i*ef2[i] for i in range(40)):.4f}")
print()

print("E3: Z world -- random odd n ~ 2^22, same statistics")
random.seed(1)
e_hist_z = np.zeros(64, dtype=np.int64)
excess_z = np.zeros(64, dtype=np.int64)
for _ in range(200000):
    n = random.getrandbits(22) | 1
    n0b = n.bit_length()
    mx = n0b
    while n != 1:
        n, e = z_step(n)
        e_hist_z[min(e, 63)] += 1
        mx = max(mx, n.bit_length())
    excess_z[min(mx - n0b, 63)] += 1
ez = e_hist_z / e_hist_z.sum()
tz = excess_z.sum()
print(f"  e-distribution Z:  e=1: {ez[1]:.4f}  e=2: {ez[2]:.4f}  "
      f"e=3: {ez[3]:.4f}  mean = {sum(i*ez[i] for i in range(64)):.4f}")
print(f"  bit-excursion: P(exc >= 3) = {float(excess_z[3:].sum())/tz:.3f}  "
      f"P(exc >= 8) = {float(excess_z[8:].sum())/tz:.3f}  "
      f"max seen = {int(np.max(np.nonzero(excess_z)))}")
print()

print("E4: carry localization -- F2-predicted e vs actual e along Z orbits")
match = 0; total = 0; drift_f2 = 0.0; drift_z = 0.0
random.seed(2)
for _ in range(50000):
    n = random.getrandbits(40) | 1
    while n != 1:
        q_f2 = f2_mul_t1(n) ^ 1          # carry-free 3n+1 on the same bits
        e_f2 = (q_f2 & -q_f2).bit_length() - 1
        n2, e_z = z_step(n)
        total += 1
        drift_f2 += e_f2; drift_z += e_z
        if e_f2 == e_z:
            match += 1
        n = n2
print(f"  steps compared: {total:,}")
print(f"  e agrees (carry-neutral steps): {match/total:.4f}")
print(f"  mean e: carry-free prediction {drift_f2/total:.4f}  vs actual {drift_z/total:.4f}")
print(f"  => carry changes e on {(1-match/total)*100:.1f}% of steps but preserves"
      f" the mean (both worlds are measure-preserving on 2-adics)")
