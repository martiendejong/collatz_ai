"""
181_n2_law.py
=============
RUNG t12: THE n^2 EXCURSION LAW, DERIVED.

Hand discovery: with I(m) = KL( Geom(1/m) || Geom(1/2) ) in bits (the LD cost
of tilting the halving-exponent mean to m), the following EXACT identity holds:

    I'(m) = log2( 2(m-1)/m )          (clean closed form)
    I(4/3) = log2(3) - 4/3            (exact equality)
    => min over m of  I(m) / (log2 3 - m)  = 1,  attained at m* = 4/3.

Consequences (theory):
 (a) cheapest climb has mean e = 4/3, i.e. P(e=1) = 1/m* = 3/4 of climb steps;
 (b) cost per climbed bit = 1  =>  P(excursion >= E bits) ~ 2^-E;
 (c) among n < 2^K expected max excursion E* ~ K  =>  peak ~ n^2:
     the DERIVATION of the empirical n^2 law of Collatz path records
     (Roosendaal: peak < ~8 n^2, conjectured sharp order n^2).

 T1: verify the identity and the minimizer numerically.
 T2: excursion record scan n < 2^20: record curve E*(K) vs K (slope -> 1),
     log2(peak)/log2(n) for record holders (-> 2).
 T3: e=1 frequency during record climbs: prediction 3/4 = 0.75
     (measured in script 179: 0.688-0.821, mean 0.765 -- retrodicted!).
"""
from math import log2

def v2(x):
    return (x & -x).bit_length() - 1

def T(n):
    q = 3 * n + 1
    e = v2(q)
    return q >> e, e

TH = log2(3)

def I(m):          # KL(Geom(1/m) || Geom(1/2)) in bits
    return log2(2 / m) + (m - 1) * log2(2 * (m - 1) / m)

print("T1: the identity")
print(f"    I(4/3)        = {I(4/3):.10f}")
print(f"    log2(3) - 4/3 = {TH - 4/3:.10f}")
gmin, mmin = 9, 0
m = 1.01
while m < TH - 1e-6:
    g = I(m) / (TH - m)
    if g < gmin:
        gmin, mmin = g, m
    m += 0.0001
print(f"    min_m I(m)/(log2 3 - m) = {gmin:.7f} at m = {mmin:.4f} "
      f"(theory: 1 at 4/3);  I'(4/3) = log2(2*(1/3)/(4/3)) = -1 exactly")

print()
print("T2: excursion records n < 2^20 (odd-iterate peaks):")
best = 0.0
records = []
LIM = 1 << 20
for n in range(3, LIM, 2):
    m, peak = n, n
    while m >= n and m != 1:
        m, _ = T(m)
        if m > peak:
            peak = m
    # full peak needs the whole orbit only if it dips below n then re-climbs
    # above old peak -- rare; do the full orbit for candidates:
    if log2(peak / n) > best - 2:
        m2, p2 = n, n
        while m2 != 1:
            m2, _ = T(m2)
            p2 = max(p2, m2)
        peak = p2
    E = log2(peak / n)
    if E > best:
        best = E
        records.append((n, peak, E))
print("    n          peak            E=log2(peak/n)  E/log2(n)  "
      "log2(peak)/log2(n)")
for n, peak, E in records[3:]:
    print(f"    {n:<10} {peak:<15} {E:10.2f}     {E/log2(n):8.3f}   "
          f"{log2(peak)/log2(n):8.3f}")

print()
print("T3: e=1 frequency during record climbs (prediction 3/4):")
for n0, _, _ in records[-6:]:
    m, peak = n0, n0
    seq = [n0]
    while m != 1:
        m, _ = T(m)
        seq.append(m)
        peak = max(peak, m)
    ip = seq.index(peak)
    m, e1, tot = n0, 0, 0
    for i in range(ip):
        m, e = T(m)
        tot += 1
        e1 += (e == 1)
    if tot:
        print(f"    n0={n0:<10} climb={tot:>3} steps  freq(e=1) = {e1/tot:.3f}")
