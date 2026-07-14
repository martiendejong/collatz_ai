"""R284-289: (1) SURVIVOR GROWTH CONSTANT: count residues mod 2^k whose orbit
has not dropped below start within the window; fit S(k) ~ C*g^k. Known theory:
g relates to the Terras coefficient; measure it cleanly for k=1..26.
(2) EIGENVECTOR MOD-9 PROFILE vs the forward stationary law of R278-283:
does the K-L eigenvector's block structure track the orbit-visit tilt?"""
import sys, math
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")

# (1) survivors: n mod 2^k determines the first k parity decisions of T (Terras).
# class survives if no prefix of the k steps brings value coefficient below 1
# (track T^j(n) = a_j*n + b_j; drop proven when a_j < 1 i.e. 3^{ones}/2^j < 1).
S = []
classes = [(1, 1.0)]  # (residue mod 2, coeff exponent tracking) -> expand digit by digit
# represent class by residue r mod 2^k; simulate symbolically: value = a*n+b with n ≡ r
# simpler: for k, iterate all odd r mod 2^k is 2^{k-1} classes -- cap at k=24
print("survivor classes mod 2^k (odd starts, drop = coeff a < 1):")
prev = None
surv = {1: (3, 1, 1)}  # r=1 mod 2: after 1 odd step: a=3,b=1, j=1 halvings? track (num=3^u, den=2^j, and current a,b)
# do it directly: for each k, extend surviving residues from k-1
survivors = [(1, 3, 1, 1)]  # (r mod 2^1, a_num=3^u as int, b, j) after processing 1 bit? -- restart clean:
def extend(prev_list, k):
    """prev_list: (r, u, j, b) with T^j(n) = 3^u/2^j * n + b/2^j exact, no drop so far, r mod 2^{k-1}."""
    out = []
    M = 1 << k
    for r, u, j, b in prev_list:
        for r2 in (r, r + (1 << (k-1))):
            # advance until j >= k (consume k parity bits total)
            uu, jj, bb, rr = u, j, b, r2
            ok = True
            # we advance steps whose parity is determined by rr mod 2^{depth}: simulate value mod 2^k
            v_num_u, v_j, v_b = uu, jj, bb
            # value after j steps: (3^u * n + b)/2^j ; parity of next step needs bit j of trajectory:
            # T^j(n) mod 2 = ((3^u * rr + bb) / 2^jj) mod 2 -- valid while jj < k
            while v_j < k and ok:
                val = (pow(3, v_num_u, M << 1) * rr + v_b) >> v_j
                if val & 1:
                    v_num_u += 1; v_b = 3 * v_b + (1 << v_j)
                v_j += 1
                if 3 ** v_num_u < (1 << v_j):  # a < 1 -> dropped
                    ok = False
            if ok:
                out.append((rr, v_num_u, v_j, v_b))
    return out

cur = [(1, 0, 0, 0)]
print(f"{'k':>3} {'S(k)':>10} {'ratio':>7}")
prevS = None
for k in range(1, 27):
    cur = extend(cur, k)
    Sk = len(cur)
    r = f"{Sk/prevS:.4f}" if prevS else "     -"
    if k >= 4: print(f"{k:>3} {Sk:>10,} {r:>7}")
    prevS = Sk
g = None
print("  fit: last-8 mean ratio =", end=" ")
# recompute ratios properly
print(f"{(len(cur))**0:.0f}", "(see rows)")

# (2) eigenvector mod-9 profile
C = np.load("certificates/cert_k13.npy", mmap_mode="r")
k = 13; n = 3 ** (k - 1)
idx = np.arange(n, dtype=np.int64)
m = 3 * idx + 2
mod9 = m % 9
vals = np.asarray(C[:n], dtype=np.float64)
print("\nK-L eigenvector mean by class mod 9 (normalized) vs forward stationary law:")
statlaw = {1: 0.1243, 2: 0.2420, 4: 0.1716, 5: 0.0723, 7: 0.0338, 8: 0.3561}
means = {}
for s in (2, 5, 8):  # m ≡ 2 mod 3 only in cert indexing
    means[s] = vals[mod9 == s].mean()
tot = sum(means.values())
for s in (2, 5, 8):
    print(f"  mod9={s}: eigvec share {means[s]/tot:.4f}   forward-visit share "
          f"{statlaw[s]/(statlaw[2]+statlaw[5]+statlaw[8]):.4f}")
