"""
180_binomial_normal_form.py
===========================
RUNG t11: the BINOMIAL NORMAL FORM of climbs (discovered by hand on the 703
wave: 32525 = 2^15 - 3^5, next 12197 = 3*2^12 - 91, then 2287 = 9*2^8 - 17).

CLAIM (exact renewal law): write n = u*2^s - v with v odd, 0 < v, and
v2(3v-1) < s (subtrahend small). Then
    T(n) = 3u * 2^{s-e} - v'   with  e = v2(3v-1),  v' = (3v-1)/2^e.
The subtrahend evolves by the 3x-1 map (= Collatz on -v), i.e. THE SHADOW
TARGET FOLLOWS THE NEGATIVE ORBIT -- this is the exact climb engine.
Special case v = 1: v' = 1 forever, e = 1 -- the classical ones-run/Mersenne
climb (t7 macro). General v: any negative integer's orbit can be tracked.
Episode ends when the 2-power budget s is exhausted (parts collide).

 V1: verify the renewal law exactly on 100k random (u, s, v).
 V2: normal-form decode of the 703 climb wave (t=15..34): best small-v target
     and its predicted successor; check v-chain = 3x-1 orbit.
 V3: is re-entry structured or Haar-luck? Along random orbits, quality
     Q = max_{v odd <= 63} [v2(n+v) - bitlen(v)]; compare tail of Q with the
     2^{-q} Haar baseline, and negative-TREE targets vs ALL small odd targets
     (prediction: per-step no difference -- tree-ness only matters for the
     PERSISTENCE of the episode, since any odd v works one step).
"""
import random
from math import log2
from collections import Counter

def v2(x):
    return (x & -x).bit_length() - 1

def T(n):
    q = 3 * n + 1
    e = v2(q)
    return q >> e, e

# ---------------- V1: renewal law ------------------------------------------
random.seed(3)
cnt = 0
for _ in range(100000):
    s = random.randint(6, 40)
    u = random.getrandbits(20) | 1
    v = (random.getrandbits(s - 3) | 1)
    if v <= 0:
        continue
    e = v2(3 * v - 1)
    if e >= s:
        continue
    n = u * (1 << s) - v
    if n <= 0:
        continue
    vp = (3 * v - 1) >> e
    nxt, ea = T(n)
    assert ea == e and nxt == 3 * u * (1 << (s - e)) - vp, (u, s, v)
    cnt += 1
print(f"V1: renewal law T(u*2^s - v) = 3u*2^(s-e) - (3v-1)/2^e verified on "
      f"{cnt:,} cases  OK")

# ---------------- V2: 703 wave in normal form --------------------------------
def best_target(n, vmax=1 << 14):
    """Best subtrahend: odd v <= vmax maximizing excess = v2(n+v) - bitlen(v)."""
    best = (-99, None, 0)
    # candidates: v = (2^d - n) mod 2^d for each depth d (the unique v per d)
    for d in range(2, n.bit_length() + 14):
        v = ((1 << d) - n) % (1 << d)
        if v == 0 or v % 2 == 0 or v > vmax:
            continue
        depth = v2(n + v)
        exc = depth - v.bit_length()
        if exc > best[0]:
            best = (exc, v, depth)
    return best

print("V2: normal-form decode of the 703 climb wave:")
n = 703
seq = [n]
while n != 1:
    n, e = T(n)
    seq.append(n)
print("    t    n          target -v   depth  excess   predicted next v")
for t in range(15, 35):
    n = seq[t]
    exc, v, depth = best_target(n)
    e = v2(3 * v - 1)
    vp = (3 * v - 1) >> e
    print(f"    {t:<4} {n:<10} -{v:<9} {depth:<6} {exc:<7} -{vp}")

# ---------------- V3: re-entry structure test --------------------------------
NEG_TREE = set()
for c in (1, 5, 7, 17, 25, 37, 55, 41, 61, 91):
    NEG_TREE.add(c)
# extend: preimages of these under the 3x-1 subtrahend map within <=63
for v in range(1, 64, 2):
    w = v
    for _ in range(30):
        w = (3 * w - 1) >> v2(3 * w - 1)
        if w in NEG_TREE:
            NEG_TREE.add(v)
            break
        if w > 1 << 20:
            break
ALL = [v for v in range(1, 64, 2)]
TREE = sorted(NEG_TREE & set(ALL))
print(f"V3: small odd targets: {len(ALL)} total, {len(TREE)} in negative tree "
      f"(3x-1 orbit stays small)")

def quality(n, targets):
    return max(v2(n + v) - v.bit_length() for v in targets)

qa = Counter()
qt = Counter()
steps = 0
random.seed(4)
for _ in range(1000):
    n = random.getrandbits(24) | 1
    while n != 1:
        qa[quality(n, ALL)] += 1
        qt[quality(n, TREE)] += 1
        steps += 1
        n, _ = T(n)
print(f"    {steps:,} orbit steps; tail P(Q >= q) vs Haar 2^-q shape:")
print("    q    all-targets    tree-only     ratio all/2^-q-fit")
ca = sum(qa.values())
tail_a = 0.0
# fit constant at q=4
Pa4 = sum(c for q, c in qa.items() if q >= 4) / ca
for q in range(4, 13):
    Pa = sum(c for qq, c in qa.items() if qq >= q) / ca
    Pt = sum(c for qq, c in qt.items() if qq >= q) / ca
    pred = Pa4 * 2 ** (-(q - 4))
    print(f"    {q:<4} {Pa:11.6f}   {Pt:11.6f}   {Pa/pred:8.3f}")
