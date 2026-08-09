# 358: testing the BALANCED (Sturmian) composition of the 186-billion ticket
# WITHOUT running the orbit — via the Stern-Brocot word algebra.
# Endpoint fixed point of a v-pattern: n = c/D, D = 2^N - 3^K, with the
# concatenation law c(UV) = c(U)*3^T(V) + 2^L(U)*c(V); single letter v:
# (T, L, c) = (1, v, 1). The balanced pattern is the mechanical word of slope
# p/q, p = N-K, q = K (v_i = 1 + sturmian bit). Mechanical words factor along
# the Stern-Brocot tree; with run-length acceleration the full algebra of the
# 72-billion-letter word costs O(CF length) operations per prime.
from math import gcd

PRIMES = [2305843009213693951, 4611686018427387847, 4611686018427387817,
          4611686018427387787, 4611686018427387761, 4611686018427387749,
          4611686018427387643, 4611686018427387619]

def make_ops(pr):
    def concat(U, V):
        Tu, Lu, cu = U
        Tv, Lv, cv = V
        return (Tu + Tv, Lu + Lv, (cu * pow(3, Tv, pr) + pow(2, Lu, pr) * cv) % pr)
    def wpow(U, e):
        R = (0, 0, 0)
        B = U
        while e:
            if e & 1:
                R = concat(R, B)
            B = concat(B, B)
            e >>= 1
        return R
    return concat, wpow

def mechanical_walk(p, q, pr):
    # exact Stern-Brocot walk to p/q (0<p<q, gcd=1) with run-length acceleration.
    # left anchor a/b < p/q < c_/d right anchor; invariant words Wl, Wr.
    # mediant word M(mediant) = concat(Wl, Wr) (validated below).
    concat, wpow = make_ops(pr)
    a, b, Wl = 0, 1, (1, 1, 1)    # letter '0': v=1
    c_, d, Wr = 1, 1, (1, 2, 1)   # letter '1': v=2
    prefixes = []
    while True:
        mp, mq = a + c_, b + d
        if mp * q == p * mq:
            return concat(Wl, Wr), prefixes
        A = p * b - a * q          # > 0
        B = c_ * q - p * d         # > 0
        if mp * q < p * mq:
            # go right: step s valid iff s < A/B, so batch t = (A-1)//B >= 1
            t = (A - 1) // B
            a, b = a + t * c_, b + t * d
            Wl = concat(Wl, wpow(Wr, t))
        else:
            # go left: step s valid iff s < B/A, so batch t = (B-1)//A >= 1
            t = (B - 1) // A
            c_, d = t * a + c_, t * b + d
            Wr = concat(wpow(Wl, t), Wr)
        prefixes.append((Wl, (a, b), Wr, (c_, d)))

def brute_algebra(p, q, pr):
    T, Lg, c = 0, 0, 0
    for i in range(q):
        v = 1 + ((i + 1) * p // q - i * p // q)
        c = (c * 3 + pow(2, Lg, pr)) % pr
        T, Lg = T + 1, Lg + v
    return T, Lg, c

# careful: brute must use same orientation as concat law; rebuild brute via concat:
def brute_algebra2(p, q, pr):
    concat, _ = make_ops(pr)
    W = (0, 0, 0)
    for i in range(q):
        v = 1 + ((i + 1) * p // q - i * p // q)
        W = concat(W, (1, v, 1))
    return W

import random
random.seed(358)
pr0 = PRIMES[0]
ok = 0
for _ in range(400):
    q = random.randrange(2, 130)
    p = random.randrange(1, q)
    if gcd(p, q) != 1:
        continue
    Wm, _ = mechanical_walk(p, q, pr0)
    Wb = brute_algebra2(p, q, pr0)
    assert Wm == Wb, (p, q, Wm, Wb)
    ok += 1
print(f"word-algebra validated against brute force on {ok} random slopes: OK")

def crt(residues, mods):
    x, M = 0, 1
    for r, m in zip(residues, mods):
        if M == 1:
            x, M = r % m, m
            continue
        t = ((r - x) * pow(M % m, -1, m)) % m
        x += M * t
        M *= m
    return x, M

def fixed_point_check(K, N, label):
    p, q = N - K, K
    g = gcd(p, q)
    res = []
    for pr in PRIMES:
        concat, wpow = make_ops(pr)
        if p == 0:
            W = wpow((1, 1, 1), q)
        elif p == q:
            W = wpow((1, 2, 1), q)
        elif g == 1:
            W, _ = mechanical_walk(p, q, pr)
        else:
            Wp, _ = mechanical_walk(p // g, q // g, pr)
            W = wpow(Wp, g)
        assert W[0] == K and W[1] == N, (W[0], W[1], K, N)
        D = (pow(2, N, pr) - pow(3, K, pr)) % pr
        res.append((W[2] * pow(D, -1, pr)) % pr)
    x, M = crt(res, PRIMES)
    xs = x if x <= M // 2 else x - M
    if abs(xs) < 2**150:
        print(f"  {label}: INTEGER fixed point n = {xs}")
        return xs
    print(f"  {label}: geen geheel vast punt onder 2^150 (CRT ~2^{abs(xs).bit_length()}) -> WEERLEGD")
    return None

print()
print("calibration on known cycles:")
fixed_point_check(1, 2, "(K=1,N=2) balanced [trivial cycle expected]")
fixed_point_check(2, 3, "(K=2,N=3) balanced [-5 expected]")
fixed_point_check(7, 11, "(K=7,N=11) balanced [-17 is NOT balanced: expect refuted]")
print()
print("THE 186-BILLION TICKET:")
K, N = 72057431991, 114208327604
fixed_point_check(K, N, "K=72,057,431,991 N=114,208,327,604 balanced")

# rotations at Stern-Brocot stage boundaries: W = P*S -> test S*P
print()
print("rotations at SB-stage boundaries:")
p, q = N - K, K
per_prime = {}
for pr in PRIMES:
    W, prefixes = mechanical_walk(p, q, pr)
    per_prime[pr] = (W, prefixes)
nstages = len(per_prime[PRIMES[0]][1])
refuted = 0
hits = []
for j in range(nstages):
    res = []
    valid = True
    for pr in PRIMES:
        concat, wpow = make_ops(pr)
        W, prefixes = per_prime[pr]
        P = prefixes[j][0]
        Tp, Lp, cp = P
        Ts, Ls = W[0] - Tp, W[1] - Lp
        if Ts <= 0 or Tp <= 0:
            valid = False
            break
        cs = ((W[2] - cp * pow(3, Ts, pr)) * pow(pow(2, Lp, pr), -1, pr)) % pr
        crot = (cs * pow(3, Tp, pr) + pow(2, Ls, pr) * cp) % pr
        D = (pow(2, N, pr) - pow(3, K, pr)) % pr
        res.append((crot * pow(D, -1, pr)) % pr)
    if not valid:
        continue
    x, M = crt(res, PRIMES)
    xs = x if x <= M // 2 else x - M
    if abs(xs) < 2**150:
        hits.append((j, xs))
    else:
        refuted += 1
print(f"  {refuted} rotaties getest en weerlegd; treffers: {hits if hits else 'geen'}")
