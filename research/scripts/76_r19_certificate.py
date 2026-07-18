# R2751-2775 WALL-2 ATTACK: the r=19 reachability certificate.
# D = 2^31 - 3^19 = 985,222,181. Bigint-bitset DP with layer recycling
# (free each dp layer once absorbed into the prefix) to stay ~4.5GB.
# 0 unreachable mod D => no 19-cycle in the critical window, by pure
# congruence DP - extending the certificate method past r=18.
import sys, math, time
sys.stdout.reconfigure(encoding="utf-8")

r = 19
j = 31
D = 2**j - 3**r
assert D == 985222181
mask = (1 << D) - 1
t0 = time.time()

def rot(x, c):
    return ((x << c) | (x >> (D - c))) & mask

dp = [0] * j
dp[0] = 1 << (pow(3, r - 1, D) % D)
for i in range(1, r):
    coef = pow(3, r - 1 - i, D)
    dp2 = [0] * j
    pref = 0
    for Jn in range(1, j):
        pref |= dp[Jn - 1]
        dp[Jn - 1] = 0                      # recycle the layer
        if pref:
            dp2[Jn] = rot(pref, (coef * pow(2, Jn, D)) % D)
    dp = dp2
    print(f"step {i}/{r-1} done  ({time.time()-t0:.0f}s)", flush=True)

out = 0
for J in range(j):
    out |= dp[J]
    dp[J] = 0
reach = bin(out).count("1")
print(f"\nr=19, D={D:,}: 0 mod D {'REACHABLE' if out & 1 else 'BLOCKED'}", flush=True)
print(f"reach fraction {reach/D:.4f}; total time {time.time()-t0:.0f}s", flush=True)
if not (out & 1):
    print("*** CERTIFICATE: no 19-cycle in the critical window (pure congruence DP) ***")
