"""R257: SUPPLY-CHAIN STRATIGRAPHY in the sigma=547 record orbit (12235060455).
R237-249 claims long paths have stratigraphy spring -> feeder -> alternator ->
repunit -> long slide. Test on the real record: along the CLIMB phase, measure
per odd value: binary trailing-ones k, alternation score of the low bits,
and compare against 200 random orbits of similar length."""
import sys, random
sys.stdout.reconfigure(encoding="utf-8")
random.seed(257)

def orbit_odds(n):
    o = [n]; m = n
    while m != 1:
        m = 3*m+1 if m % 2 else m//2
        if m % 2 and m != 1: o.append(m)
    return o

def tail_ones(n):
    k = 0
    while n & 1: n >>= 1; k += 1
    return k

def alt_score(n, bits=12):
    # fraction of adjacent bit-pairs that alternate in the low `bits` bits
    s = bin(n)[2:][-bits:]
    if len(s) < 2: return 0.0
    return sum(1 for a, b in zip(s, s[1:]) if a != b) / (len(s)-1)

REC = 12235060455
ro = orbit_odds(REC)
peak_i = max(range(len(ro)), key=lambda i: ro[i])
climb = ro[:peak_i+1]
print(f"record orbit: {len(ro)} odd steps, peak at step {peak_i} value ~2^{ro[peak_i].bit_length()}")

k_climb = [tail_ones(v) for v in climb]
a_climb = [alt_score(v) for v in climb]
print(f"CLIMB ({len(climb)} odds): mean k = {sum(k_climb)/len(k_climb):.3f}, "
      f"mean alt-score = {sum(a_climb)/len(a_climb):.4f}, max k = {max(k_climb)}")
desc = ro[peak_i+1:]
k_d = [tail_ones(v) for v in desc]; a_d = [alt_score(v) for v in desc]
print(f"DESCENT ({len(desc)} odds): mean k = {sum(k_d)/len(k_d):.3f}, "
      f"mean alt-score = {sum(a_d)/len(a_d):.4f}")

# baseline: random orbits, same statistic over their climb phases
ks, asc = [], []
for _ in range(200):
    n = random.randrange(3, 1 << 34) | 1
    o = orbit_odds(n)
    pi = max(range(len(o)), key=lambda i: o[i])
    c = o[:pi+1]
    if len(c) < 3: continue
    ks.append(sum(tail_ones(v) for v in c)/len(c))
    asc.append(sum(alt_score(v) for v in c)/len(c))
print(f"\nbaseline random climbs: mean k = {sum(ks)/len(ks):.3f}, "
      f"mean alt-score = {sum(asc)/len(asc):.4f}")

# where in the climb do the big-k bursts sit? show the 10 largest-k moments
top = sorted(range(len(climb)), key=lambda i: -k_climb[i])[:10]
print("\nbiggest fuel moments in the record climb (step, k, alt-score of low 12 bits):")
for i in sorted(top):
    print(f"  step {i:>3}: k={k_climb[i]:>2}  alt={a_climb[i]:.3f}  low bits ...{bin(climb[i])[-16:]}")
