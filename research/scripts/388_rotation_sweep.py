# 388: ROTATION SWEEP of ticket 1 — key insight: the fixed point of the
# rotated word is the next orbit element of the rational fixed point
# n1 = c/(2^N - 3^K). So sweeping rotations = iterating
#   n_{i+1} = (3 n_i + 1)/2^{v_i}  (mod primes),  v_i = mechanical word letter.
# Any rotation with an integer fixed point (any cycle in this pattern class)
# shows up as a small CRT value. Calibration: (K,N)=(2,3) must detect -5/-7.
import time

PRIMES = [2305843009213693951, 4611686018427387847, 4611686018427387817,
          4611686018427387787, 4611686018427387761, 4611686018427387749,
          4611686018427387643, 4611686018427387619]

def word_c(K, N, pr):
    # (T, L, c) of the mechanical word via SB walk (as in 385)
    def concat(U, V):
        Tu, Lu, cu = U; Tv, Lv, cv = V
        return (Tu+Tv, Lu+Lv, (cu*pow(3, Tv, pr) + pow(2, Lu, pr)*cv) % pr)
    def wpow(U, e):
        R = (0, 0, 0); B = U
        while e:
            if e & 1: R = concat(R, B)
            B = concat(B, B); e >>= 1
        return R
    p, q = N-K, K
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

def crt_small(res, mods, bound):
    x, M = 0, 1
    for r, m in zip(res, mods):
        if M == 1: x, M = r % m, m; continue
        t = ((r-x)*pow(M % m, -1, m)) % m
        x += M*t; M *= m
    xs = x if x <= M//2 else x - M
    return xs if abs(xs) < bound else None

def sweep(K, N, R, bound, label):
    p, q = N-K, K
    n = []
    inv2 = []
    for pr in PRIMES:
        W = word_c(K, N, pr)
        assert W[0] == K and W[1] == N
        D = (pow(2, N, pr) - pow(3, K, pr)) % pr
        n.append((W[2]*pow(D, -1, pr)) % pr)
        inv2.append([pow(pow(2, v, pr), -1, pr) for v in range(4)])
    t0 = time.time()
    hits = []
    fp = 0
    for i in range(1, R+1):
        x = crt_small(n, PRIMES, bound)
        if x is not None:
            hits.append((i-1, x))
        v = 1 + ((i*p)//q - ((i-1)*p)//q)
        for j, pr in enumerate(PRIMES):
            n[j] = ((3*n[j] + 1)*inv2[j][v]) % pr
    dt = time.time() - t0
    print(f"{label}: {R:,} rotaties in {dt:.0f}s; integer-treffers: "
          f"{hits[:6] if hits else 'GEEN'}", flush=True)
    return hits

# calibration: (K,N) = (2,3) — the -5 pattern
sweep(2, 3, 2, 1 << 40, "kalibratie (K=2,N=3)")
# ticket 1: first million rotations
sweep(72057431991, 114208327604, 1000000, 1 << 160, "ticket 1 (K=72e9)")
