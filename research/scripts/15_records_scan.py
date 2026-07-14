"""E25: vectorized stopping-time records scan. Usage: python 15_records_scan.py START END
Tracks records of sigma(n) = T-steps until value < n, for n = 3 mod 4 in [START, END).
Overflow-safe: values capped; escapees re-run in bigint python."""
import sys, numpy as np
sys.stdout.reconfigure(encoding="utf-8")

START = int(sys.argv[1]) if len(sys.argv) > 1 else 3
END = int(sys.argv[2]) if len(sys.argv) > 2 else 200_000_000
RECMIN = int(sys.argv[3]) if len(sys.argv) > 3 else 100   # only report records above this
CHUNK = 1 << 22
CAP = 700
LIM = np.int64(1) << np.int64(61)

def stop_bigint(n, cap=2000):
    m = n; t = 0
    while t < cap:
        m = (3 * m + 1) // 2 if m % 2 else m // 2
        t += 1
        if m < n: return t
    return -1

best = RECMIN
records = []
s = START + ((3 - START) % 4)
while s < END:
    e = min(s + CHUNK * 4, END)
    n0 = np.arange(s, e, 4, dtype=np.int64)
    m = n0.copy()
    alive = np.ones(n0.shape, dtype=bool)
    t = 0
    while alive.any() and t < CAP:
        mm = m[alive]
        odd = (mm & 1).astype(bool)
        mm[odd] = 3 * mm[odd] + 1
        mm >>= 1
        m[alive] = mm
        t += 1
        idx = np.where(alive)[0]
        dropped = mm < n0[idx]
        over = mm > LIM
        # overflow escapees: resolve in python bigint
        for i in idx[over]:
            st = stop_bigint(int(n0[i]))
            if st > best:
                best = st; records.append((int(n0[i]), st))
                print(f"RECORD {n0[i]} sigma={st}", flush=True)
        alive[idx[dropped | over]] = False
        # record check for those dropping exactly now with t > best
        if t > best:
            for i in idx[dropped]:
                if t > best:
                    best = t; records.append((int(n0[i]), t))
                    print(f"RECORD {n0[i]} sigma={t}", flush=True)
    # survivors past CAP: bigint
    for i in np.where(alive)[0]:
        st = stop_bigint(int(n0[i]))
        if st > best:
            best = st; records.append((int(n0[i]), st))
            print(f"RECORD {n0[i]} sigma={st}", flush=True)
    s = e
print("DONE", START, END, "records:", records)
