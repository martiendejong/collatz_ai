# 348: reverse-engineering the producers of giant ones-runs.
# Forward macro-step (course 2.2 form): odd n0 = a*2^k - 1 (a odd, exact run k),
#   A = a*3^k, v = v2(A-1), m1 = (A-1)/2^v  (next odd number).
# Backward: given a target cylinder  m1 = r mod 2^M  and a chosen branch (k, v):
#   a = (r*2^v + 1) * 3^(-k) mod 2^(M+v),  n0 = a*2^k - 1 mod 2^(M+v+k).
# One residue class per branch. Questions:
#  (Q1) the classical claim: the predecessor of the all-ones number 2^K - 1
#       (branch k=1, v=1) has front a = (4^m - 1)/3 = binary 0101...01 (alternating).
#  (Q2) what leads to the alternating pattern? one more pullback: does the
#       required tail become periodic with period 6 = ord(2 mod 9), then 18, ...?
#  (Q3) full chains: verify forward that a member of the pulled-back class
#       really produces the giant ones-run after t steps.
from math import log2

def v2(x):
    c = 0
    while x % 2 == 0:
        x //= 2
        c += 1
    return c

def ones_run(x):
    c = 0
    while x % 2 == 1:
        c += 1
        x //= 2
    return c

def macro(n):  # one combined step, returns (next_odd, k, v)
    k = ones_run(n)
    a = (n + 1) >> k
    A = a * 3**k
    v = v2(A - 1)
    return (A - 1) >> v, k, v

def pullback(r, M, k, v):
    # solve for the class of n0 mod 2^(M+v+k) that maps onto m1 = r mod 2^M via branch (k,v)
    Mp = M + v
    inv3k = pow(3**k, -1, 2**Mp)
    a = ((r * 2**v + 1) * inv3k) % 2**Mp
    assert a % 2 == 1, "front a must be odd for an exact run of k"
    n0 = (a * 2**k - 1) % 2**(Mp + k)
    return n0, Mp + k, a

def bits(x, m):
    return bin(x % 2**m)[2:].zfill(m)

K = 9
target = 2**K - 1          # all-ones number: 111111111
M = K + 1                  # pin it mod 2^(K+1) so the run is exactly K
print(f"TARGET: m = {target} = {bits(target, M)} (run of {K} ones), pinned mod 2^{M}")
print()

# Q1: one step back, branch (1,1)
n1, M1, a1 = pullback(target, M, 1, 1)
print(f"t=1 back (k=1,v=1): n = {bits(n1, M1)} mod 2^{M1}   front a = {bits(a1, M1-1)}")
m4 = (4**5 - 1) // 3       # (4^m-1)/3 with 2m-1 = K+1... check the classical form
print(f"   classical (4^m-1)/3 alternating check: a mod 2^10 = {bits(a1, 10)}  vs 0101010101")
print()

# Q2/Q3: chain further back, always branch (1,1).
# Closed form: the backward map of branch (1,1) is m -> (4m-1)/3 in Z2, with fixed
# point 1, so the t-step producer class is  x_t = 1 + (4/3)^t (m - 1)  (denominator 3^t).
r, Mc = target, M
chain = []
for t in range(1, 5):
    r, Mc, a = pullback(r, Mc, 1, 1)
    chain.append((t, r, Mc, a))
    s = bits(r, Mc)
    x_t = (1 + 4**t * (target - 1) * pow(3**t, -1, 2**Mc)) % 2**Mc
    print(f"t={t} back: class mod 2^{Mc}: ...{s[-40:] if len(s)>40 else s}   "
          f"closed form 1+(4/3)^t(m-1): {'OK' if x_t == r else 'FAIL'}")
print()

# Q3: forward verification — take the smallest member of each class and run it forward
print("forward verification of the pulled-back classes:")
for t, r, Mc, a in chain:
    n = r if r % 2 == 1 else r + 2**Mc
    orig = n
    ks = []
    for _ in range(t):
        n, k, v = macro(n)
        ks.append((k, v))
    run = ones_run(n)
    ok = "OK" if run >= K else "FAIL"
    print(f"  t={t}: start {orig} -> after {t} macro steps: run of {run} ones "
          f"(target >= {K}) {ok}   branches used: {ks}")
print()

# the emblematic concrete chain, small numbers
print("emblematic chain (smallest instances):")
n = 169  # 10101001
trail = [n]
for _ in range(3):
    n, k, v = macro(n)
    trail.append(n)
print("  " + " -> ".join(f"{x}={bin(x)[2:]}" for x in trail[:3]))
