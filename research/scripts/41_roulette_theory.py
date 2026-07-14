"""R309-315: EXACT ROULETTE THEORY & THE BRIDGE TRIANGLE.
Under the geometric-w model, the forward mod-3^j chain is exactly computable:
from state r (mod 3^j, r not div 3): base b=(3r+1) mod 3^j (depends only on
r mod 3^{j-1}), next = b * 2^{-w} mod 3^j with P(w)=2^{-w}.
Closed form at j=2: pi = (8,16,11,4,2,22)/63 for classes (1,2,4,5,7,8).
Compare THEORY vs ORBIT-MEASURED vs EIGENVECTOR block means at j=2,3,4."""
import sys, random, math
import numpy as np
from collections import Counter
sys.stdout.reconfigure(encoding="utf-8")
random.seed(309)

def theory_stationary(j):
    M = 3 ** j
    states = [r for r in range(M) if r % 3 != 0]
    idx = {r: i for i, r in enumerate(states)}
    n = len(states)
    P = np.zeros((n, n))
    inv2 = pow((M + 1) // 2, 1, M)  # 2^{-1} mod 3^j
    for r in states:
        b = (3 * r + 1) % M
        p = 0.5; x = b
        for w in range(1, 80):
            x = (x * inv2) % M
            P[idx[r], idx[x]] += p
            p *= 0.5
    # renormalize rows (truncation) and get stationary
    P /= P.sum(1, keepdims=True)
    v = np.ones(n) / n
    for _ in range(2000): v = v @ P
    return states, v

# orbit measurement mod 81
stat = Counter()
for _ in range(80_000):
    n0 = random.randrange(3, 1 << 44) | 1
    m = n0
    while m != 1:
        m = 3*m+1 if m % 2 else m//2
        if m % 2 and m != 1: stat[m % 81] += 1

# eigenvector block means mod 81 (classes == 2 mod 3 only)
C = np.load("certificates/cert_k13.npy", mmap_mode="r")
k = 13; N = 3 ** (k - 1)
mvals = 3 * np.arange(N, dtype=np.int64) + 2
ev = np.asarray(C[:N], dtype=np.float64)

print("closed-form check at j=2: predicted (8,16,11,4,2,22)/63:")
states, v = theory_stationary(2)
pred = {1: 8/63, 2: 16/63, 4: 11/63, 5: 4/63, 7: 2/63, 8: 22/63}
for s in states:
    print(f"  mod9={s}: numeric {v[states.index(s)]:.4f}  closed-form {pred[s]:.4f}")

for j in (2, 3, 4):
    M = 3 ** j
    states, v = theory_stationary(j)
    th = {s: v[i] for i, s in enumerate(states)}
    # restrict to classes == 2 mod 3 (eigenvector coverage), renormalize all three
    cls = [s for s in states if s % 3 == 2]
    t = np.array([th[s] for s in cls]); t /= t.sum()
    o = np.array([stat[s % M] if M <= 81 else 0 for s in cls], dtype=float)
    # orbit counts at mod M: aggregate stat (mod 81) down to mod M
    o = np.array([sum(c for r, c in stat.items() if r % M == s) for s in cls], dtype=float)
    o /= o.sum()
    e = np.array([ev[(mvals % M) == s].mean() for s in cls]); e /= e.sum()
    r_to = np.corrcoef(t, o)[0, 1]
    r_te = np.corrcoef(t, e)[0, 1]
    r_oe = np.corrcoef(o, e)[0, 1]
    print(f"\nmod 3^{j} = {M} ({len(cls)} classes ==2 mod 3):")
    print(f"  corr(theory, orbits) = {r_to:.4f} | corr(theory, eigvec) = {r_te:.4f} | corr(orbits, eigvec) = {r_oe:.4f}")
    if j == 2:
        print("  shares (theory | orbit | eigvec):",
              " ".join(f"{s}:({t[i]:.3f}|{o[i]:.3f}|{e[i]:.3f})" for i, s in enumerate(cls)))
