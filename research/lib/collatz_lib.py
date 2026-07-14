"""Shared machinery for the Collatz family/pair research program.

Coordinates: for odd n, write n + 1 = a * 2^k with a odd.
  a = family parameter, k = sequence index (= count of trailing binary 1s of n).
The segment: k consecutive steps n -> (3n+1)/2 (each exact), ending at a*3^k - 1.
Macro-step: segment + the cascade of w = v2(a*3^k - 1) halvings to the next odd.
"""

def v2(x):
    v = 0
    while x % 2 == 0:
        x //= 2
        v += 1
    return v


def trailing_ones(n):
    k = 0
    while n & 1:
        n >>= 1
        k += 1
    return k


def coords(n):
    """odd n -> (a, k) with n + 1 = a * 2^k, a odd."""
    k = trailing_ones(n)
    return (n + 1) >> k, k


def starter(a, k):
    return a * (1 << k) - 1


def macro_step(a, k):
    """(a,k) -> (a', k', w): one full segment + halving cascade."""
    x = a * 3 ** k - 1
    w = v2(x)
    m = x >> w
    k2 = trailing_ones(m)
    return (m + 1) >> k2, k2, w


def syracuse(n):
    """collapsed odd->odd step; returns (next_odd, halvings)."""
    x = 3 * n + 1
    w = v2(x)
    return x >> w, w


def syracuse_c(n, c):
    x = c * n + 1
    w = v2(x)
    return x >> w, w


def collatz_map(n):
    return 3 * n + 1 if n % 2 else n // 2


def map_3nm1(n):
    return 3 * n - 1 if n % 2 else n // 2


def map_5np1(n):
    return 5 * n + 1 if n % 2 else n // 2


def ternary(n):
    if n == 0:
        return "0"
    d = []
    while n:
        d.append(str(n % 3))
        n //= 3
    return "".join(reversed(d))


def base_repr(n, b):
    if n == 0:
        return "0"
    d = []
    while n:
        d.append(str(n % b))
        n //= b
    return "".join(reversed(d))


def ternlen(n):
    L = 0
    while n:
        n //= 3
        L += 1
    return L


def cycle_min(n, f):
    """minimum element of the cycle eventually reached under map f (Floyd)."""
    slow = fast = n
    while True:
        slow = f(slow)
        fast = f(f(fast))
        if slow == fast:
            break
    cyc = [slow]
    m = f(slow)
    while m != slow:
        cyc.append(m)
        m = f(m)
    return min(cyc)
