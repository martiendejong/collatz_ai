# 353: the run-automaton lookup table — "which pattern above the run leads to
# how many ones at the tail next step".
# State before a conversion: tail run of k ones, front a (odd). The step is the
# SAME affine rule every time: A = a*3^k, v = v2(A-1), m1 = (A-1)/2^v.
# The next ones-count k' is DETERMINED by (k, a mod 2^w) for a window w:
# a lookup table, time-invariant — the "changes the same way every time".
# We tabulate: (k, a mod 32) -> (v, k'), marking entries undetermined within
# the window, and count how thin the classes for k' >= 3 are.
def v2(x):
    c = 0
    while x % 2 == 0:
        x //= 2
        c += 1
    return c

W = 7  # window: a mod 2^7
print(f"run-automaton table: (k, a mod 2^{W}) -> (v, k')   [window {W} bits]")
print("(k' shown as '>=x' when the window cannot determine the full run)")
for k in range(1, 4):
    rows = []
    rich = []
    for a in range(1, 2**W, 2):
        A = a * 3**k
        Am = (A - 1) % 2**(W + k)   # bits of A-1 known from a mod 2^W (approx window)
        v = v2(A - 1)
        # count next ones within the reliable window
        kp = 0
        det = True
        while (A >> (v + kp)) & 1:
            kp += 1
            if v + kp >= W - 1:   # ran out of certain bits
                det = False
                break
        rows.append((a, v, kp, det))
        if kp >= 3:
            rich.append(a)
    shown = ", ".join(f"a={a}:(v={v},k'={'>=' if not d else ''}{kp})" for a, v, kp, d in rows[:8])
    print(f"  k={k}: {shown} ...")
    print(f"       classes with k'>=3 within window: {len(rich)}/{2**(W-1)} -> {rich}")
print()
# neutrality: expected k' over uniform a, per k (should be ~2 for every k)
import random
random.seed(353)
for k in range(1, 6):
    tot = 0
    n = 40000
    for _ in range(n):
        a = random.getrandbits(48) * 2 + 1
        A = a * 3**k
        v = v2(A - 1)
        m1 = (A - 1) >> v
        kp = 0
        while m1 % 2 == 1:
            kp += 1
            m1 //= 2
        tot += kp
    print(f"  E[k' | k={k}] = {tot/n:.3f}   (neutral value 2)")
