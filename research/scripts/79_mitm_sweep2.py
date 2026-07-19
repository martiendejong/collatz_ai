# R3641+: MITM certificates S=36..41. Per-shape enumeration keeps memory
# bounded (max one binomial bucket at a time) for c >= 31.
import sys, math, time, bisect
from array import array
sys.stdout.reconfigure(encoding="utf-8")

def enum_shape(L, s_target, M, offset_pow):
    """array of W mod M over all words of length L with exactly s_target rises.
    Positions processed descending; rise at i adds 3^(rises_above)*2^(i+offset).
    offset_pow: precomputed pow2[i] table (len L) already including offset."""
    pow3 = [pow(3, s, M) for s in range(s_target+1)]
    out = array('q')
    # iterative DFS with pruning: remaining positions must fit remaining rises
    stack = [(L-1, 0, 0)]
    while stack:
        i, s, W = stack.pop()
        if s == s_target:
            out.append(W)
            continue
        if i < 0:
            continue
        need = s_target - s
        if need > i + 1:
            continue
        stack.append((i-1, s, W))
        stack.append((i-1, s+1, (W + pow3[s]*offset_pow[i]) % M))
    return out

def certify(S, D):
    M = 2**D - 3**S
    c = D // 2
    L2 = D - c
    t0 = time.time()
    pow2_pre = [pow(2, i, M) for i in range(c)]
    pow2_suf = [pow(2, c + i, M) for i in range(L2)]  # suffix carries 2^c folded in
    hits = 0
    for s in range(max(0, S - L2), min(S, c) + 1):
        s2 = S - s
        T1 = enum_shape(c, s, M, pow2_pre)
        T1 = array('q', sorted(T1))
        # need: 3^s2*W1 + W2' == 0 with W2' = 2^c*W2 (folded): W1 == -W2'*inv(3^s2)
        inv = pow(pow(3, s2, M), -1, M)
        # enumerate suffix shape s2, look up
        pow3 = [pow(3, t, M) for t in range(s2+1)]
        stack = [(L2-1, 0, 0)]
        while stack:
            i, t, W = stack.pop()
            if t == s2:
                need = (-W * inv) % M
                j = bisect.bisect_left(T1, need)
                if j < len(T1) and T1[j] == need: hits += 1
                continue
            if i < 0: continue
            if s2 - t > i + 1: continue
            stack.append((i-1, t, W))
            stack.append((i-1, t+1, (W + pow3[t]*pow2_suf[i]) % M))
        del T1
    return hits, time.time() - t0, M

if __name__ == "__main__":
    for S in (36, 37, 38, 39, 40, 41):
        D = math.ceil(S * math.log2(3))
        h, t, M = certify(S, D)
        print(f"S={S:3d} D={D:3d} M={M}  cycle-words={h}  [{t:.0f}s]", flush=True)
