"""
150_audit_prior_claims.py
==========================
AUDIT of the auxiliary Diophantine claims from earlier sessions, before
building on them. Claims under audit:

 (1) "PPD Obstruction Theorem": if a prime p | 2^B - 3^A has large order
     (ord_p(2/3) > A), then no Collatz cycle with parameters (A,B) exists.
     -> Suspicion: VACUOUS/INVALID. In the exact integer identity
        n1*(2^B - 3^A) = C, every prime dividing 2^B - 3^A divides C
        automatically. p | C is not an extra condition; an obstruction
        would need to show v_p(C) < v_p(2^B - 3^A) for ALL admissible
        d-sequences, which was never established.

 (2) "Primitive prime divisors of 2^B - 3^A have ord_p(2/3) = B."
     -> Testable directly by factoring.

 (3) "Wieferich-like pairs": pairs (A,B) where every prime p | 2^B-3^A has
     ord_p(3/2) <= A. Claimed to be 'the only pairs surviving PPD'.
     -> Test how often the condition holds at small scale, and whether
        rad(2^B-3^A) | G := gcd(3^A-2^A, 2^(B-A)-1) characterizes them.

 (4) "n1 < 5x10^103 from Baker" -- textual audit (see conclusions below).

Standard-form conventions used here and in scripts 151-152:
  cycle of V(n) = (3n+1)/2^e on odd n, with A odd steps, B = sum(e_i);
  n1*(2^B - 3^A) = C,  C = sum_{i=1}^{A} 3^(A-i) * 2^(E_{i-1}),
  0 = E_0 < E_1 < ... < E_{A-1} <= B-1.
"""
import math, random

# ---------- factoring (Pollard rho) ----------
def is_prime(n):
    if n < 2: return False
    for p in (2,3,5,7,11,13,17,19,23,29,31,37):
        if n % p == 0: return n == p
    d, s = n-1, 0
    while d % 2 == 0: d //= 2; s += 1
    for a in (2,3,5,7,11,13,17,19,23,29,31,37):
        x = pow(a, d, n)
        if x in (1, n-1): continue
        for _ in range(s-1):
            x = x*x % n
            if x == n-1: break
        else:
            return False
    return True

def pollard(n):
    if n % 2 == 0: return 2
    while True:
        x = random.randrange(2, n); y = x; c = random.randrange(1, n); d = 1
        while d == 1:
            x = (x*x + c) % n
            y = (y*y + c) % n; y = (y*y + c) % n
            d = math.gcd(abs(x-y), n)
        if d != n: return d

def factor(n):
    if n == 1: return {}
    if is_prime(n): return {n: 1}
    d = pollard(n)
    f = factor(d)
    for p, e in factor(n // d).items():
        f[p] = f.get(p, 0) + e
    return f

def mult_order(a, p):
    """ord of a mod prime p (a coprime to p)."""
    o = p - 1
    for q in factor(p - 1):
        while o % q == 0 and pow(a, o // q, p) == 1:
            o //= q
    return o

THETA = math.log2(3)

print("="*78)
print("AUDIT (1)+(2)+(3): factor D = 2^B - 3^A for near-convergent pairs, A<=34")
print("="*78)
rows = []
wieferich_like = []
prim_ord_eq_B = 0; prim_total = 0
for A in range(2, 35):
    B = math.ceil(A * THETA)
    D = (1 << B) - 3**A
    if D <= 0:
        B += 1; D = (1 << B) - 3**A
    fac = factor(D)
    all_small_order = True
    parts = []
    for p in sorted(fac):
        if p < 5:
            parts.append(f"{p}^{fac[p]}" if fac[p] > 1 else f"{p}")
            continue
        inv2 = pow(2, -1, p)
        o32 = mult_order(3 * inv2 % p, p)     # ord_p(3/2)
        o23 = mult_order(2 * pow(3, -1, p) % p, p)  # ord_p(2/3) (same thing)
        parts.append(f"{p}(ord32={o32})")
        prim_total += 1
        if o32 == B: prim_ord_eq_B += 1
        if o32 > A: all_small_order = False
    G = math.gcd(3**A - 2**A, (1 << (B - A)) - 1)
    radD_divides_G = all((p < 3) or (G % p == 0) for p in fac)
    tag = ""
    if all_small_order:
        wieferich_like.append((A, B))
        tag = "  <-- Wieferich-like (all ord_p(3/2) <= A)"
    print(f"A={A:>2} B={B:>2} D={D:<16} G={G:<8} radD|G={str(radD_divides_G):<5} "
          f"{' '.join(parts)}{tag}")

print()
print(f"pairs with ALL prime factors of ord_p(3/2) <= A ('Wieferich-like'): "
      f"{wieferich_like if wieferich_like else 'NONE in range'}")
print(f"primes p|D with ord_p(3/2) == B exactly: {prim_ord_eq_B}/{prim_total} "
      f"(claim (2) said 'always B' for primitive primes)")

print()
print("="*78)
print("AUDIT (1): is p | C an 'obstruction'?  Direct check.")
print("="*78)
print("""
In the exact identity n1*(2^B - 3^A) = C, take ANY prime p | 2^B - 3^A.
Then p divides the left-hand side, hence p | C -- AUTOMATICALLY, with no
condition on ord_p(anything). The 'PPD Obstruction Theorem' as recorded
('p | 2^B-3^A with ord_p(2/3) > A implies no cycle') asserts an obstruction
where none was derived: it never showed that v_p(C) < v_p(2^B-3^A) for all
admissible E-sequences. VERDICT: claim (1) is NOT a theorem as stated.
""")

# demonstrate: for a sample pair, MANY admissible sequences have C ≡ 0 mod p
A, B = 12, 20
D = (1 << B) - 3**A
p = max(factor(D))
from itertools import combinations
hits = tot = 0
for E in combinations(range(1, B), A - 1):
    C = pow(3, A-1) + sum(3**(A-1-i) * (1 << E[i-1]) for i in range(1, A))
    tot += 1
    if C % p == 0: hits += 1
print(f"demo (A={A},B={B}): largest prime of D is p={p}; among all {tot} admissible")
print(f"E-sequences, {hits} satisfy C == 0 (mod p)  (expected ~ tot/p = {tot//p})")
print("=> single-prime divisibility of C is abundant, never an obstruction per se.")

print()
print("="*78)
print("AUDIT (3)+(4): conditional and textual verdicts")
print("="*78)
print("""
(3) The k>=3 statement ('if rad(D) | G and D | G^k then k >= 3 for A >= 100')
    is arithmetically sound as a CONDITIONAL: G <= 2^(B-A)-1 < 2^(0.586A),
    so G^2 < 2^(1.18A) < D ~ 2^(1.585A - o(A)) by Baker/Rhin. But it
    constrains only the class 'rad(D) | G', which the table above shows is
    empty at small scale and was never proven to contain all cycle-admissible
    pairs (that inference rested on the invalid claim (1)).
    VERDICT: conditionally true, currently without force.

(4) 'n1 < 5x10^103': no absolute upper bound on cycle elements follows from
    Baker-type results, because A is unbounded. Bounds of that shape exist
    only for cycles with a BOUNDED NUMBER OF BLOCKS (m-cycles, Simons-de
    Weger). VERDICT: unfounded as stated; discard.
""")
