"""R296-301: CONTROLS. The doorway machinery must transfer to 3n-1 (structure)
while the CONCLUSION differs (three cycles). Measure for 3n-1:
(1) basin shares of cycles {1}, {5}, {17} among random odd starts;
(2) the entry doorway of each cycle: distribution of the odd value at which the
    orbit first enters the cycle, and the caste (mod 3) quantization of the last
    pre-cycle odd -> verifies our laws are MAP-STRUCTURAL, not 3n+1-magic;
(3) 5n+1 divergence check with gate analog: does the (2^k-1)/5 'gate family'
    exist and is it caste-blocked the same way? (2^k == 1 mod 5 iff 4|k)."""
import sys, random
from collections import Counter
sys.stdout.reconfigure(encoding="utf-8")
random.seed(296)


# build cycle sets properly by iteration
def cyc_of(seed):
    s, m = set(), seed
    while m not in s:
        s.add(m)
        m = 3*m-1 if m % 2 else m//2
    return frozenset(x for x in s if x % 2)
C1, C5, C17 = cyc_of(1), cyc_of(5), cyc_of(17)
print("3n-1 odd cycles:", sorted(C1), sorted(C5), sorted(C17))

N = 100_000
basin = Counter(); door = Counter()
for _ in range(N):
    n = random.randrange(3, 1 << 40) | 1
    m = n; prev_odd = None; last_odd = n
    seen_guard = 0
    while True:
        if m % 2 and m in C1: c = 1; break
        if m % 2 and m in C5: c = 5; break
        if m % 2 and m in C17: c = 17; break
        nm = 3*m-1 if m % 2 else m//2
        if nm % 2: prev_odd, last_odd = last_odd, nm
        m = nm
    basin[c] += 1
    door[(c, prev_odd)] += 1

tot = sum(basin.values())
print(f"basin shares: cycle1 {basin[1]/tot:.4f}, cycle5 {basin[5]/tot:.4f}, cycle17 {basin[17]/tot:.4f}")
print("\ntop doors per cycle (penultimate odd before entering):")
for c in (1, 5, 17):
    ds = Counter({v: n for (cc, v), n in door.items() if cc == c})
    top = ds.most_common(4)
    tt = sum(ds.values())
    print(f"  cycle {c}: " + ", ".join(f"{v} ({n/tt:.3f}, mod3={v%3})" for v, n in top if v))

# (3) 5n+1 gate family: gates g with 5g+1 = 2^k -> g=(2^k-1)/5, integer iff k%4==0
print("\n5n+1 gate family (2^k-1)/5, k=4,8,12,...:", [(2**k-1)//5 for k in (4, 8, 12, 16)])
print("  caste mod 5 of gates:", [((2**k-1)//5) % 5 for k in (4, 8, 12, 16)],
      "-> multiples of 5 are the sourceless caste in 5n+1;")
g5 = [(2**k-1)//5 for k in (4, 8, 12, 16, 20, 24)]
print("  which gates are multiples of 5 (shut doors):", [g % 5 == 0 for g in g5])
