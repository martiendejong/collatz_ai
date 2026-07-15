"""R542-550: THE HOP TAX, exact. An orbit shadows phase -1 at depth k
(= k trailing ones). After the full burn, the refill k' is set by the clockwork.
The set {n : first r successive refills all >= j} is a finite union of residue
classes mod 2^K -- count it EXACTLY (theorem-grade, like the census).
Prediction from memorylessness: density = (2^{-(j-1)})^r."""
import sys
sys.stdout.reconfigure(encoding="utf-8")

def refill_chain(n, r):
    """follow r macro-burns, return list of refill depths (trailing ones)."""
    out = []
    for _ in range(r):
        k = 0
        while (n >> k) & 1: k += 1
        out.append(k)
        a = (n + 1) >> k
        v = a * 3**k - 1
        while v % 2 == 0: v //= 2
        n = v
    return out

K = 22  # count over odd residues mod 2^K exactly (chain depths stay well below K bits of influence... 
# NOTE: refills depend on bits beyond 2^K for some n; use direct integer computation on all odd n < 2^K:
M = 1 << K
print(f"exact count over all odd n < 2^{K}:")
print(f"{'j':>2} {'r':>2} {'density':>10} {'pred (2^-(j-1))^r':>17} {'ratio':>6}")
import itertools
# efficient: single pass computing chains of depth 4
counts = {}
for j in (2, 3, 4):
    for r in (1, 2, 3):
        counts[(j, r)] = 0
tot = 0
for n in range(3, M, 2):
    ch = refill_chain(n, 3)
    tot += 1
    for j in (2, 3, 4):
        ok = True
        for r in (1, 2, 3):
            if ch[r-1] < j: ok = False
            if ok: counts[(j, r)] += 1
for j in (2, 3, 4):
    for r in (1, 2, 3):
        d = counts[(j, r)]/tot
        p = (2.0**(-(j-1)))**r
        print(f"{j:>2} {r:>2} {d:>10.6f} {p:>17.6f} {d/p:>6.3f}")
