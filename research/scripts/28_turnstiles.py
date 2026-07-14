"""R254: THE TURNSTILES. Claim from R252-253: almost all orbits funnel through
penultimate odd 13 or 53 (the two open low rungs of gate 5's odd-w ladder;
rungs w=1 (->3) and w=7 (->213) are springs, blocked). Verify directly and
map the funnel one layer deeper (antepenultimate odd)."""
import sys, random
from collections import Counter
sys.stdout.reconfigure(encoding="utf-8")
random.seed(254)

N = 200_000
pen, ante = Counter(), Counter()
for _ in range(N):
    n = random.randrange(3, 1 << 42) | 1
    odds = []
    m = n
    while m != 1:
        m = 3*m+1 if m % 2 else m//2
        if m % 2 and m != 1: odds.append(m)
    o = [n] + odds
    pen[o[-2] if len(o) >= 2 else None] += 1
    ante[o[-3] if len(o) >= 3 else None] += 1

print("penultimate odd (the turnstile):")
for v, c in pen.most_common(8):
    print(f"  {v}: {c/N:.4f}   mod3={v%3 if v else '-'}")
print(f"  cumulative top-2: {sum(c for _,c in pen.most_common(2))/N:.4f}")
print("\nantepenultimate odd (one layer up the funnel):")
for v, c in ante.most_common(12):
    print(f"  {v}: {c/N:.4f}   mod3={v%3 if v else '-'}")
