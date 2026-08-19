# 398: the bias channel of the halving law. The diagonal correlation
# corr(p_t, n_t) = bias of the junk f_t; it halves per two steps. Hypothesis:
# each halving comes from ONE specific fresh AND entering the bias channel.
# Test: condition the bias on products n_a & n_b; if one product restores the
# bias by exactly 2x, the channel is a chain of specific ANDs (proof target
# for the carry recursion). Light sampling (memory-friendly).
import numpy as np
rng = np.random.default_rng(398)
M = 400000
n0 = (rng.integers(0, 2**62, M, dtype=np.int64) << 1) | 1
B = 18
bits = np.array([((n0 >> b) & 1) for b in range(B)], dtype=np.int8)
n = n0.copy()
print(f"{'t':>3} {'bias(p^n_t)':>11} {'beste conditionering (a&b)':>26} {'bias|ab=1':>9} {'bias|ab=0':>9}")
for t in range(1, 13):
    odd = (n & 1).astype(bool)
    n = np.where(odd, 3*n + 1, n) >> 1
    if t < 4 or t >= B:
        continue
    res = ((n & 1).astype(np.int8) ^ bits[t])   # p_t XOR n_t = f_t
    s = 1 - 2*res.astype(np.float64)
    bias = abs(float(s.mean()))
    best = (0, 0, 0.0, 0.0, 0.0)
    for a in range(max(0, t-6), t):
        for b2 in range(a+1, t+1):
            ab = (bits[a] & bits[b2]).astype(bool)
            if ab.sum() < 1000: continue
            b1v = abs(float(s[ab].mean()))
            b0v = abs(float(s[~ab].mean()))
            gain = max(b1v, b0v)
            if gain > best[4]:
                best = (a, b2, b1v, b0v, gain)
    a, b2, b1v, b0v, _ = best
    print(f"{t:>3} {bias:>11.4f} {'n%d & n%d' % (a, b2):>26} {b1v:>9.4f} {b0v:>9.4f}", flush=True)
