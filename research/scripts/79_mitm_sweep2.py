# R3641+: MITM certificates S=35..41. Per-shape enumeration bounds memory
# to one binomial bucket; numpy handles sorting and batched lookups
# (8 bytes/element - the sorted()-of-Python-ints MemoryError is gone).
import sys, math, time
import numpy as np
from array import array
sys.stdout.reconfigure(encoding="utf-8")

def enum_shape(L, s_target, M, pow2):
    """array('q') of W mod M over length-L words with exactly s_target rises.
    Positions descending; rise at position i adds 3^(rises_above)*pow2[i]."""
    pow3 = [pow(3, s, M) for s in range(s_target+1)]
    out = array('q')
    stack = [(L-1, 0, 0)]
    while stack:
        i, s, W = stack.pop()
        if s == s_target:
            out.append(W)
            continue
        if i < 0 or s_target - s > i + 1:
            continue
        stack.append((i-1, s, W))
        stack.append((i-1, s+1, (W + pow3[s]*pow2[i]) % M))
    return out

def certify(S, D, chunk=4_000_000):
    M = 2**D - 3**S
    c = D // 2
    L2 = D - c
    t0 = time.time()
    pow2_pre = [pow(2, i, M) for i in range(c)]
    pow2_suf = [pow(2, c + i, M) for i in range(L2)]  # 2^c folded into suffix
    hits = 0
    for s in range(max(0, S - L2), min(S, c) + 1):
        s2 = S - s
        T1raw = enum_shape(c, s, M, pow2_pre)
        T1 = np.frombuffer(T1raw, dtype=np.int64).copy()
        del T1raw
        T1.sort()
        inv = pow(pow(3, s2, M), -1, M)
        pow3 = [pow(3, t, M) for t in range(s2+1)]
        buf = array('q')
        def flush():
            nonlocal hits
            if not len(buf): return
            a = np.frombuffer(buf, dtype=np.int64)
            idx = np.searchsorted(T1, a)
            ok = idx < len(T1)
            hits += int(np.sum(T1[idx[ok]] == a[ok]))
        stack = [(L2-1, 0, 0)]
        while stack:
            i, t, W = stack.pop()
            if t == s2:
                buf.append((-W * inv) % M)
                if len(buf) >= chunk:
                    flush(); buf = array('q')
                continue
            if i < 0 or s2 - t > i + 1:
                continue
            stack.append((i-1, t, W))
            stack.append((i-1, t+1, (W + pow3[t]*pow2_suf[i]) % M))
        flush()
        del T1, buf
    return hits, time.time() - t0, M

if __name__ == "__main__":
    for S in (35, 36, 37, 38, 39, 40, 41):
        D = math.ceil(S * math.log2(3))
        h, t, M = certify(S, D)
        print(f"S={S:3d} D={D:3d} M={M}  cycle-words={h}  [{t:.0f}s]", flush=True)
