"""R278-283: CASTE MARKOV DYNAMICS along forward odd orbits.
Theory: next odd = (3m+1)/2^w with 3m+1 == 1 mod 3, so next mod 3 = 2^w mod 3 =
1 if w even, 2 if w odd -> caste = w-parity, and if w ~ geometric(1/2) indep of
caste: P(next=2)=2/3, P(next=1)=1/3, NO memory. Test mod 3, then refine mod 9
(6 states 1,2,4,5,7,8): does mod-9 have genuine Markov memory? Compare stationary
distribution with uniform and measure the spectral gap of the empirical chain."""
import sys, random
import numpy as np
from collections import Counter
sys.stdout.reconfigure(encoding="utf-8")
random.seed(278)

trans3 = Counter(); trans9 = Counter(); stat9 = Counter()
N = 60_000
for _ in range(N):
    n = random.randrange(3, 1 << 44) | 1
    m = n; prev3 = prev9 = None
    while m != 1:
        nm = 3*m+1 if m % 2 else m//2
        if nm % 2 and nm != 1:
            c3, c9 = nm % 3, nm % 9
            if prev3 is not None:
                trans3[(prev3, c3)] += 1
                trans9[(prev9, c9)] += 1
            stat9[c9] += 1
            prev3, prev9 = c3, c9
        m = nm

print("mod-3 caste transition matrix P(next | cur):  [prediction: rows equal (1/3, 2/3)]")
for a in (1, 2):
    tot = sum(trans3[(a, b)] for b in (1, 2))
    print(f"  from {a}: " + "  ".join(f"P({b})={trans3[(a,b)]/tot:.4f}" for b in (1, 2)))

S = [1, 2, 4, 5, 7, 8]
print("\nmod-9 stationary distribution [uniform = 0.1667]:")
tot = sum(stat9.values())
print("  " + "  ".join(f"{s}:{stat9[s]/tot:.4f}" for s in S))

P = np.zeros((6, 6))
for i, a in enumerate(S):
    row = sum(trans9[(a, b)] for b in S)
    for j, b in enumerate(S):
        P[i, j] = trans9[(a, b)]/row if row else 0
print("\nmod-9 transition matrix (rows=from):")
for i, a in enumerate(S):
    print(f"  {a}: " + " ".join(f"{P[i,j]:.3f}" for j in range(6)))
ev = np.sort(np.abs(np.linalg.eigvals(P)))[::-1]
print(f"\n|eigenvalues|: {[f'{e:.4f}' for e in ev]}")
print(f"spectral gap 1-|lambda_2| = {1-ev[1]:.4f}  (memory half-life = {np.log(0.5)/np.log(ev[1]):.2f} steps)"
      if ev[1] > 0 else "no memory")
# chi^2 test row-vs-row: is there ANY mod-9 memory?
rows_equal = np.allclose(P, P[0], atol=0.01)
print(f"rows near-identical (no memory) at 0.01: {rows_equal}")
