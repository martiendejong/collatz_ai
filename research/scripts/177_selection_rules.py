"""
177_selection_rules.py
======================
Machine verification of Obs 339-340 (rigidity bridge t8-t9) + first measurements.

 V1: exit law (l=1 branch): K' = v2(3^K a + 1) - 1 and a' = (3^K a + 1)/2^{K'+1}
 V2: Mersenne dead-stop: a = 1, K odd  =>  K' = 1  (LTE)
 V3: channel table mod 8 (c = 3^K a mod 8):
     c=3: l=1, K'=1 forced   c=7: l=1, K'>=2   c=5: l=2   c=1: l>=3
 V4: decode record orbits (27, 703, 26623, 626331, 837799) as macro/channel
     chains; verify EVERY deep restart (K'>=2 after l=1 exit) sits in channel 7,
     and print the 27-chain against the hand-decoding of Obs 339.
 E5: excursion-shadow test (Obs 337 prediction): e=1 frequency during the climb
     phase of record orbits; prediction: > 0.586 while climbing.
 E6: WHERE DOES THE ARITHMETIC HIDE? Joint law of (l, K') under Haar-random
     input: test whether (l, K') are independent Geom(1/2) x Geom(1/2).
     (If yes, the symbol-level law is Haar-blind and the killer must live in
     the deterministic evolution of a's residues -- the automaton of t9.)
"""
import random
from math import log2

def v2(x):
    return (x & -x).bit_length() - 1

def T(n):
    q = 3 * n + 1
    e = v2(q)
    return q >> e, e

def macro(n):
    """n odd -> (K, a, c, l, u): run length K, multiplier a, channel c,
    exit valuation l, landing odd u. Macro consumes K odd steps of T."""
    K = v2(n + 1)
    a = (n + 1) >> K
    x = 3**K * a
    l = v2(x - 1)
    u = (x - 1) >> l
    c = x & 7
    return K, a, c, l, u

# ---------------- V0: macro consistency against direct iteration ----------
random.seed(7)
for _ in range(20000):
    n = random.getrandbits(40) | 1
    K, a, c, l, u = macro(n)
    m = n
    for _ in range(K):
        m, e = T(m)
    assert m == u, (n, u, m)
print("V0: macro(n) == T^K(n) on 20,000 random 40-bit odd n  OK")

# ---------------- V1: exit law for l=1 ------------------------------------
cnt = 0
for _ in range(200000):
    K = v2(random.getrandbits(20) | 1 << 19) + 1   # rough geometric K >= 1
    a = random.getrandbits(40) | 1
    x = 3**K * a
    if v2(x - 1) != 1:
        continue
    u = (x - 1) >> 1
    Kp = v2(u + 1)
    assert Kp == v2(x + 1) - 1
    assert (u + 1) >> Kp == (x + 1) >> (Kp + 1)
    cnt += 1
print(f"V1: exit law K' = v2(3^K a + 1) - 1 verified on {cnt:,} l=1 cases  OK")

# ---------------- V2: Mersenne dead-stop -----------------------------------
for K in range(1, 200, 2):
    x = 3**K
    assert v2(x - 1) == 1            # l = 1
    u = (x - 1) >> 1
    assert v2(u + 1) == 1, K         # K' = 1
print("V2: Mersenne dead-stop (a=1, K odd => l=1, K'=1) verified K=1..199  OK")

# ---------------- V3: channel table ---------------------------------------
seen = {1: [], 3: [], 5: [], 7: []}
for _ in range(200000):
    K = v2(random.getrandbits(20) | 1 << 19) + 1
    a = random.getrandbits(40) | 1
    x = 3**K * a
    c = x & 7
    l = v2(x - 1)
    u = (x - 1) >> l
    Kp = v2(u + 1)
    seen[c].append((l, Kp))
ok = all(l == 1 and kp == 1 for l, kp in seen[3])
ok &= all(l == 1 and kp >= 2 for l, kp in seen[7])
ok &= all(l == 2 for l, kp in seen[5])
ok &= all(l >= 3 for l, kp in seen[1])
print(f"V3: channel table (c=3: l=1,K'=1 | c=7: l=1,K'>=2 | c=5: l=2 | c=1: l>=3)"
      f"  {'OK' if ok else 'FAILED'}")

# ---------------- V4: record orbit decoding --------------------------------
def decode(n0, verbose=False):
    n = n0
    chain = []
    while n != 1:
        K, a, c, l, u = macro(n)
        Kp = v2(u + 1) if u != 1 else 0
        deep = (l == 1 and Kp >= 2)
        # selection rule: deep <=> channel 7
        assert deep == (c == 7) or u == 1, (n, K, a, c, l, u)
        chain.append((n, K, a, c, l, u))
        if verbose:
            grow = K * log2(3) - (K + l)
            tag = "DEEP" if deep else ("dead" if c == 3 else "")
            print(f"    n={n:<8} K={K:<2} a={a:<7} c={c} l={l:<2} -> u={u:<8}"
                  f"  {'+' if grow>0 else ''}{grow:.2f} bits  {tag}")
        n = u
    return chain

print("V4: record orbits as macro/channel chains (selection rule asserted "
      "at every step):")
print("  n0 = 27 (hand-decoding of Obs 339):")
decode(27, verbose=True)
for n0 in (703, 26623, 626331, 837799):
    ch = decode(n0)
    deep = sum(1 for (_, K, a, c, l, u) in ch if c == 7)
    print(f"  n0={n0}: {len(ch)} macro-steps, {deep} channel-7 (deep) hits, "
          f"selection rule OK on all")

# ---------------- E5: excursion-shadow test --------------------------------
print("E5: e=1 frequency during climb phase (prediction: > 0.586 to climb):")
for n0 in (27, 703, 26623, 626331, 837799):
    # find peak over odd iterates
    n, peak = n0, n0
    seq = [n0]
    while n != 1:
        n, e = T(n)
        seq.append(n)
        peak = max(peak, n)
    # climb = steps until first time peak is reached
    n, e1, tot = n0, 0, 0
    while n != peak:
        n, e = T(n)
        tot += 1
        e1 += (e == 1)
    print(f"  n0={n0:<7} peak={peak:<10} climb={tot:>3} odd-steps  "
          f"freq(e=1)={e1/tot:.3f}  net={(log2(3)*tot - sum(0 for _ in []) ):.0f}")

# ---------------- E6: joint law (l, K') under Haar --------------------------
print("E6: joint law of (l, K') under Haar input (independence test):")
from collections import Counter
joint = Counter()
N = 500000
for _ in range(N):
    a = random.getrandbits(48) | 1
    K = v2(random.getrandbits(20) | 1 << 19) + 1
    x = 3**K * a
    l = v2(x - 1)
    u = (x - 1) >> l
    Kp = v2(u + 1)
    joint[(min(l, 5), min(Kp, 5))] += 1
print("      P(l,K') / [P(l)P(K')]  (1.00 = independent):")
Pl = Counter(); Pk = Counter()
for (l, k), c in joint.items():
    Pl[l] += c; Pk[k] += c
hdr = "      l\\K' " + "".join(f"{k:>7}" for k in range(1, 6))
print(hdr)
for l in range(1, 6):
    row = f"      {l:<5}"
    for k in range(1, 6):
        p = joint.get((l, k), 0) / N
        pi = (Pl[l] / N) * (Pk[k] / N)
        row += f"{p/pi:7.3f}" if pi > 0 else "      -"
    print(row)
