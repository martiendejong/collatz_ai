"""R1491-1505: THE CONGRUENCE OBSTRUCTION HUNT. Census (script 72) saw 0
nontrivial hits vs ~9.7 expected (P ~ 6e-5): W mod D does NOT equidistribute.
Hypothesis: for each r there is a prime q | D with W != 0 mod q for ALL words
-> a finite congruence CERTIFICATE excluding r-cycles in the critical window.
DP: achievable set of W mod q over all compositions (j_1..j_r) of j."""
import sys, math
import numpy as np
import sympy
sys.stdout.reconfigure(encoding="utf-8")
ALPHA = math.log2(3)

def achievable_zero(r, j, q):
    """can W = sum_{i=0}^{r-1} 3^{r-1-i} 2^{J_i} be 0 mod q, J_0=0<J_1<...<J_{r-1}<=j-1?"""
    if q == 1: return True
    # dp[J] = boolean array over w mod q, after adding terms 0..i with J_i = J
    dp = [np.zeros(q, dtype=bool) for _ in range(j)]
    dp[0][pow(3, r - 1, q) % q] = True          # i=0 term, J_0 = 0
    for i in range(1, r):
        coef = pow(3, r - 1 - i, q)
        ndp = [np.zeros(q, dtype=bool) for _ in range(j)]
        # prefix OR over previous J
        pref = np.zeros(q, dtype=bool)
        for Jp in range(j - 1):
            pref |= dp[Jp]
            Jn = Jp + 1
            if Jn > j - 1: break
            if pref.any():
                c = (coef * pow(2, Jn, q)) % q
                ndp[Jn] |= np.roll(pref, c)
        # allow J_n > Jp+1: pref already accumulates all smaller J; but roll depends on Jn:
        # redo properly: for each Jn, OR over all Jp < Jn of dp[Jp], then roll by coef*2^Jn
        dp2 = [np.zeros(q, dtype=bool) for _ in range(j)]
        pref = np.zeros(q, dtype=bool)
        for Jn in range(1, j):
            pref |= dp[Jn - 1]
            if pref.any():
                c = (coef * pow(2, Jn, q)) % q
                dp2[Jn] = np.roll(pref, c)
        dp = dp2
    out = np.zeros(q, dtype=bool)
    for J in range(j): out |= dp[J]
    return bool(out[0]), int(out.sum())

print(f"{'r':>3s} {'j':>3s} {'D factorization':40s} obstruction?")
found_all = True
for r in range(3, 25):
    j = math.ceil(r * ALPHA)
    if 2 ** j - 3 ** r <= 0: j += 1
    D = 2 ** j - 3 ** r
    fac = sympy.factorint(D)
    fs = "*".join(f"{p}^{e}" if e > 1 else str(p) for p, e in fac.items())
    blocked = None
    detail = []
    for p in sorted(fac):
        if p > 200000: detail.append(f"{p}:skip"); continue
        ok0, nach = achievable_zero(r, j, p)
        detail.append(f"{p}:{'0-REACHABLE' if ok0 else 'BLOCKED'}({nach}/{p})")
        if not ok0 and blocked is None:
            blocked = p
    tag = f"YES q={blocked}" if blocked else "no small-factor block"
    if blocked is None: found_all = False
    print(f"{r:3d} {j:3d} {fs:40s} {tag}   [{'; '.join(detail)}]")
print(f"\nall r blocked by a small prime factor: {found_all}")
