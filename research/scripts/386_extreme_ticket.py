# 386: THE MOST EXTREME TESTABLE TICKET — push the balanced-pattern test to
# cycle candidates of astronomical length. Cost of the Stern-Brocot word
# algebra is O(log K), so the testability frontier is set only by the
# precision of log2(3) and the CRT prime count. Targets: tickets of length
# ~1e20, ~1e100, ~1e300 steps (the last: vastly more steps than atoms in
# the observable universe ~1e80).
from mpmath import mp, mpf, log
from math import gcd
import time

mp.dps = 800
ALPHA = log(3)/log(2)
LN2 = log(2); LN3 = log(3)

# continued fraction of alpha; collect convergents
x = ALPHA
p0, q0, p1, q1 = 1, 0, int(ALPHA), 1
convs = [(p1, q1)]
frac = x - int(x)
while q1 < 10**310:
    x = 1/frac
    a = int(x)
    frac = x - a
    p0, p1 = p1, a*p1 + p0
    q0, q1 = q1, a*q1 + q0
    convs.append((p1, q1))

def is_prime(n):
    for p in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]:
        if n % p == 0:
            return n == p
    d, s = n-1, 0
    while d % 2 == 0: d //= 2; s += 1
    for a in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]:
        y = pow(a, d, n)
        if y in (1, n-1): continue
        for _ in range(s-1):
            y = y*y % n
            if y == n-1: break
        else:
            return False
    return True

def primes61(count):
    out = []
    c = (1 << 61) - 1
    while len(out) < count:
        if is_prime(c): out.append(c)
        c -= 2
    return out

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

def walk(p, q, pr):
    concat, wpow = make_ops(pr)
    a, b, Wl = 0, 1, (1, 1, 1)
    c_, d, Wr = 1, 1, (1, 2, 1)
    while True:
        if (a+c_)*q == p*(b+d):
            return concat(Wl, Wr)
        A = p*b - a*q; B = c_*q - p*d
        if (a+c_)*q < p*(b+d):
            t = (A-1)//B
            a, b = a+t*c_, b+t*d
            Wl = concat(Wl, wpow(Wr, t))
        else:
            t = (B-1)//A
            c_, d = t*a+c_, t*b+d
            Wr = concat(wpow(Wl, t), Wr)

def crt(residues, mods):
    x, M = 0, 1
    for r, m in zip(residues, mods):
        if M == 1: x, M = r % m, m; continue
        t = ((r-x)*pow(M % m, -1, m)) % m
        x += M*t; M *= m
    return x, M

for target in [10**20, 10**100, 10**300]:
    # first upper-side convergent with q >= target
    tick = None
    for p_, q_ in convs:
        if q_ >= target and mpf(p_)/q_ > ALPHA:
            tick = (q_, p_)
            break
    if tick is None:
        continue
    K, N = tick
    Delta = N*LN2 - K*LN3
    Vstar = K/(3*Delta)
    l2V = float(log(Vstar)/LN2)
    nprimes = int((2*l2V + 128)/61) + 1
    t0 = time.time()
    PR = primes61(nprimes)
    p_sl, q_sl = N-K, K
    assert gcd(p_sl, q_sl) == 1
    res = []
    for pr in PR:
        W = walk(p_sl, q_sl, pr)
        assert W[0] == K and W[1] == N
        D = (pow(2, N, pr) - pow(3, K, pr)) % pr
        res.append((W[2]*pow(D, -1, pr)) % pr)
    x, M = crt(res, PR)
    xs = x if x <= M//2 else x - M
    bound = 1 << (int(2*l2V) + 64)
    verdict = "GEHEEL VAST PUNT!" if abs(xs) < bound else "WEERLEGD"
    dt = time.time() - t0
    print(f"ticket lengte ~1e{len(str(K+N))-1}: K = 1e{len(str(K))-1}-schaal "
          f"({str(K)[:12]}...), lengte K+N = {str(K+N)[:12]}... "
          f"({len(str(K+N))} cijfers)")
    print(f"  venster Delta ~ 1e{int(mp.log10(Delta))}; vereiste elementen "
          f"~2^{l2V:.0f}; {nprimes} priemgetallen")
    print(f"  gebalanceerd patroon: {verdict} (CRT ~2^{abs(xs).bit_length()}, "
          f"grens 2^{int(2*l2V)+64}) in {dt:.1f}s", flush=True)
    print()
