# R3591+: MITM certificates for S=33..35, array-based (memory-lean).
import sys, math, time, bisect
from array import array
sys.stdout.reconfigure(encoding="utf-8")

def enum_halves_arr(L, M):
    res = {}
    pow2 = [pow(2, i, M) for i in range(L)]
    pow3 = [pow(3, s, M) for s in range(L+1)]
    stack = [(L-1, 0, 0)]
    while stack:
        i, s, W = stack.pop()
        if i < 0:
            res.setdefault(s, array('q')).append(W)
            continue
        stack.append((i-1, s, W))
        stack.append((i-1, s+1, (W + pow3[s]*pow2[i]) % M))
    for s in res:
        res[s] = array('q', sorted(res[s]))
    return res

def certify(S, D):
    M = 2**D - 3**S
    c = D // 2
    t0 = time.time()
    H1 = enum_halves_arr(c, M)
    H2 = enum_halves_arr(D - c, M)
    hits = 0
    for s in range(0, S+1):
        s2 = S - s
        if s not in H1 or s2 not in H2: continue
        f = (-pow(2, c, M) * pow(pow(3, s2, M), -1, M)) % M
        T1 = H1[s]
        for W2 in H2[s2]:
            need = (f * W2) % M
            j = bisect.bisect_left(T1, need)
            if j < len(T1) and T1[j] == need: hits += 1
    return hits, time.time() - t0, M

if __name__ == "__main__":
    for S in (33, 34, 35):
        D = math.ceil(S * math.log2(3))
        h, t, M = certify(S, D)
        print(f"S={S:3d} D={D:3d} M={M}  cycle-words={h}  [{t:.0f}s]", flush=True)
