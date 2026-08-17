# 388b: background — 100M rotations of ticket 1 (continuation of 388)
import time, sys
sys.path.insert(0, r"E:\projects\collatz\research\scripts")
exec(open(r"E:\projects\collatz\research\scripts\388_rotation_sweep.py").read().split("# calibration")[0])

K, N = 72057431991, 114208327604
p, q = N-K, K
n = []; inv2 = []
for pr in PRIMES:
    W = word_c(K, N, pr)
    D = (pow(2, N, pr) - pow(3, K, pr)) % pr
    n.append((W[2]*pow(D, -1, pr)) % pr)
    inv2.append([pow(pow(2, v, pr), -1, pr) for v in range(4)])
R = 100_000_000
bound = 1 << 160
t0 = time.time()
hits = []
CHK = 10_000_000
for i in range(1, R+1):
    x = crt_small(n, PRIMES, bound)
    if x is not None:
        hits.append((i-1, x))
        print("HIT:", i-1, x, flush=True)
    v = 1 + ((i*p)//q - ((i-1)*p)//q)
    for j, pr in enumerate(PRIMES):
        n[j] = ((3*n[j] + 1)*inv2[j][v]) % pr
    if i % CHK == 0:
        print(f"{i:,} rotaties, {time.time()-t0:.0f}s, treffers: {len(hits)}", flush=True)
print(f"KLAAR: {R:,} rotaties in {time.time()-t0:.0f}s; treffers: {hits if hits else 'GEEN'}")
