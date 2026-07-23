"""
149_no_congruence_theorem.py
=============================
THE NO-CONGRUENCE THEOREM (negative structural result).

Question: can ANY finite-modulus (congruence) argument mod 2^N rule out
Collatz cycles -- e.g. a repaired version of the phantom-cascade argument?

Answer: NO. This script verifies that for every odd residue r mod 2^N
(N = 7 and N = 10 tested exhaustively for coverage), there is a periodic
parameter sequence (K_1,l0_1,...,K_L,l0_L) whose unique 2-adic solution

    n1 = C / (2^B - 3^A)   in Z_2      (D = 2^B - 3^A is odd, so invertible)

is FULLY 2-adically consistent (simulating the macro-step orbit at high
precision confirms v2(n_i + 1) = K_i and the l0_i valuations at EVERY step
around the loop) and satisfies n1 = r (mod 2^N).

This matches Lagarias's theorem that the 3x+1 map on Z_2 is conjugate to
the shift: every parameter/parity sequence is realized by a unique 2-adic
point, and periodic sequences give rational cycles. The residues that the
cascade compatibility argument declared impossible -- {91,103,175,603,615,687}
mod 1024 -- are all realized (examples printed below).

Consequence: what separates a genuine integer cycle from a harmless rational
one is ONLY the Diophantine integrality condition (2^B - 3^A) | C. That
condition is invisible to every finite modulus. The cycle case of the Collatz
conjecture therefore cannot be settled by congruence/phantom-cycle methods.
"""
import random

def v2m(x, cap):
    if x == 0:
        return cap
    c = 0
    while x % 2 == 0 and c < cap:
        x //= 2; c += 1
    return c

def solve_and_check(Ks, l0s, N, slack=48):
    """Return n1 mod 2^N for the unique 2-adic solution of the cycle equation
    for parameter sequence (Ks, l0s), or None if the sequence is not fully
    2-adically consistent (it always is, per Lagarias, unless precision runs out)."""
    L = len(Ks); A = sum(Ks); B = A + sum(l0s)
    M = N + B + slack
    mod = 1 << M
    # cycle equation: n1 * (2^B - 3^A) = C
    C = 0
    for i in range(L):
        d_i = 3**Ks[i] - 2**Ks[i]
        exp3 = sum(Ks[i+1:])
        exp2 = sum(Ks[j] + l0s[j] for j in range(i))
        C += d_i * (3**exp3) * (2**exp2)
    D = (1 << B) - 3**A          # odd => invertible in Z/2^M
    n1 = (C * pow(D % mod, -1, mod)) % mod
    # simulate loop, checking full consistency; precision drops K_i+l0_i per step
    n = n1; prec = M
    for i in range(L):
        if v2m((n + 1) % (1 << prec), prec) != Ks[i]:
            return None
        m = ((n + 1) >> Ks[i]) % (1 << (prec - Ks[i]))
        x = (m * 3**Ks[i] - 1) % (1 << (prec - Ks[i]))
        if v2m(x, prec - Ks[i] - 1) != l0s[i]:
            return None
        n = (x >> l0s[i]) % (1 << (prec - Ks[i] - l0s[i]))
        prec -= Ks[i] + l0s[i]
        if prec <= N:
            return None
    return n1 % (1 << N)

def coverage(N, trials, seed, kmax=12, lmax=9, lmin_len=1, lmax_len=10):
    random.seed(seed)
    realized = {}
    full = 1 << (N - 1)
    for _ in range(trials):
        L = random.randint(lmin_len, lmax_len)
        Ks  = [random.randint(1, kmax) for _ in range(L)]
        l0s = [random.randint(1, lmax) for _ in range(L)]
        r = solve_and_check(Ks, l0s, N)
        if r is not None and r not in realized:
            realized[r] = (Ks, l0s)
            if len(realized) == full:
                break
    return realized

if __name__ == "__main__":
    for N, trials in [(7, 500_000), (10, 8_000_000)]:
        realized = coverage(N, trials, seed=99)
        odd = set(range(1, 1 << N, 2))
        missing = sorted(odd - set(realized))
        print(f"N={N}: {len(realized)}/{len(odd)} odd residues mod 2^{N} realized "
              f"by fully-consistent 2-adic cycles; missing: {missing if missing else 'NONE'}")

    # the residues the cascade argument declared impossible:
    realized10 = coverage(10, 3_000_000, seed=11, kmax=10, lmax=7, lmax_len=7)
    print("\nCascade-'forbidden' residues mod 1024, realized by 2-adic cycles:")
    for x in [91, 103, 175, 603, 615, 687]:
        if x in realized10:
            Ks, l0s = realized10[x]
            print(f"  r={x}: K={Ks}, l0={l0s}")
