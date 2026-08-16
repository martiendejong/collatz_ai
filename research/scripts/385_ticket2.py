# 385: the NEXT cycle candidate — ticket 2: K = 137,528,045,312 (the t=2
# mediant 6,586,818,670 + 2*65,470,613,321), N = 217,976,794,617, length
# 355,504,839,929. Exact window/price + balanced-pattern refutation via the
# Stern-Brocot word algebra (same machinery as script 358 / Obs 545).
from math import gcd
from decimal import Decimal, getcontext
getcontext().prec = 80

K = 137528045312
N = 217976794617
ln2 = Decimal(2).ln(); ln3 = Decimal(3).ln()
Delta = N*ln2 - K*ln3
Vstar = K/(3*Delta)
print(f"ticket 2: K = {K:,}  N = {N:,}  lengte = {K+N:,}")
print(f"venster Delta = {Delta:.4E}; vereiste elementschaal V* = {Vstar:.3E} "
      f"= 2^{float((Vstar.ln()/ln2)):.2f}")
print(f"mediant-check: 6586818670 + 2*65470613321 = {6586818670 + 2*65470613321:,}")
print()

PRIMES = [2305843009213693951, 4611686018427387847, 4611686018427387817,
          4611686018427387787, 4611686018427387761, 4611686018427387749,
          4611686018427387643, 4611686018427387619]

def make_ops(pr):
    def concat(U, V):
        Tu, Lu, cu = U; Tv, Lv, cv = V
        return (Tu+Tv, Lu+Lv, (cu*pow(3, Tv, pr) + pow(2, Lu, pr)*cv) % pr)
    def wpow(U, e):
        R = (0, 0, 0); B = U
        while e:
            if e & 1: R = concat(R, B)
            B = concat(B, B); e >>= 1
        return R
    return concat, wpow

def mechanical_walk(p, q, pr):
    concat, wpow = make_ops(pr)
    a, b, Wl = 0, 1, (1, 1, 1)
    c_, d, Wr = 1, 1, (1, 2, 1)
    prefixes = []
    while True:
        mp, mq = a+c_, b+d
        if mp*q == p*mq:
            return concat(Wl, Wr), prefixes
        A = p*b - a*q; B = c_*q - p*d
        if mp*q < p*mq:
            t = (A-1)//B
            a, b = a+t*c_, b+t*d
            Wl = concat(Wl, wpow(Wr, t))
        else:
            t = (B-1)//A
            c_, d = t*a+c_, t*b+d
            Wr = concat(wpow(Wl, t), Wr)
        prefixes.append((Wl, (a, b), Wr, (c_, d)))

def crt(residues, mods):
    x, M = 0, 1
    for r, m in zip(residues, mods):
        if M == 1:
            x, M = r % m, m; continue
        t = ((r-x)*pow(M % m, -1, m)) % m
        x += M*t; M *= m
    return x, M

p, q = N-K, K
g = gcd(p, q)
print(f"helling p/q = {p:,}/{q:,}  gcd = {g}")
res = []
for pr in PRIMES:
    concat, wpow = make_ops(pr)
    if g == 1:
        W, _ = mechanical_walk(p, q, pr)
    else:
        Wp, _ = mechanical_walk(p//g, q//g, pr)
        W = wpow(Wp, g)
    assert W[0] == K and W[1] == N
    D = (pow(2, N, pr) - pow(3, K, pr)) % pr
    res.append((W[2]*pow(D, -1, pr)) % pr)
x, M = crt(res, PRIMES)
xs = x if x <= M//2 else x - M
if abs(xs) < 2**160:
    print(f"gebalanceerd patroon: GEHEEL vast punt n = {xs}  <-- !!!")
else:
    print(f"gebalanceerd patroon: geen geheel vast punt onder 2^160 "
          f"(CRT-waarde ~2^{abs(xs).bit_length()}) -> WEERLEGD")
# rotations at SB-stage boundaries
per_prime = {}
for pr in PRIMES:
    W, prefixes = mechanical_walk(p//g, q//g, pr) if g > 1 else mechanical_walk(p, q, pr)
    per_prime[pr] = (W, prefixes)
nst = len(per_prime[PRIMES[0]][1])
refuted = hits = 0
for j in range(nst):
    resr = []; ok = True
    for pr in PRIMES:
        concat, wpow = make_ops(pr)
        W, prefixes = per_prime[pr]
        if g > 1:
            W = make_ops(pr)[1](W, g)
        P = prefixes[j][0]
        Tp, Lp, cp = P
        Ts, Ls = W[0]-Tp, W[1]-Lp
        if Ts <= 0 or Tp <= 0:
            ok = False; break
        cs = ((W[2] - cp*pow(3, Ts, pr))*pow(pow(2, Lp, pr), -1, pr)) % pr
        crot = (cs*pow(3, Tp, pr) + pow(2, Ls, pr)*cp) % pr
        D = (pow(2, N, pr) - pow(3, K, pr)) % pr
        resr.append((crot*pow(D, -1, pr)) % pr)
    if not ok: continue
    x, M = crt(resr, PRIMES)
    xs = x if x <= M//2 else x - M
    if abs(xs) < 2**160: hits += 1
    else: refuted += 1
print(f"rotaties op SB-grenzen: {refuted} weerlegd, {hits} treffers")
