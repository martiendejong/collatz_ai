"""Obs 618: phase words as family coordinates; alternators are phase-extremal.

(1) Exact bijection (provable, verified): a mod 2^j <-> (coset sign,
    phase r_j) via a*3^r = +-1 mod 2^j -- so phase words are EXACTLY
    equidistributed over families; they are coordinates, not statistics.
(2) Alternator extremality: A_y = (4^y-1)/3 satisfies A_y*3 + 1 = 4^y,
    so v2(A_y*3^k + 1) at k=1 equals 2y: the deepest possible fuel
    event sits at the FIRST clock tick; phase r_j = 1 for all j <= 2y.
    The gateways of 1 are the maximally fuel-resonant families.
"""

def v2(x):
    v = 0
    while x % 2 == 0:
        x //= 2
        v += 1
    return v

def phase(a, j):
    """r in [0, 2^(j-2)) with a*3^r = +-1 mod 2^j (sign by coset)."""
    sign = 1 if a % 8 in (1, 3) else -1
    mod = 1 << j
    tgt = (sign) % mod
    r, p = 0, 1  # p = 3^r mod 2^j
    per = 1 << (j - 2)
    for r in range(per):
        if (a * p) % mod == tgt:
            return r
        p = (3 * p) % mod
    raise AssertionError(f"no phase for a={a}, j={j}")

if __name__ == "__main__":
    # (1) bijection: for each j, the map a -> r_j is exactly 4-to-1 per
    # residue (2 cosets x the a mod 2^j lift freedom collapses); check
    # exact equidistribution of r_j over odd a mod 2^j.
    from collections import Counter
    for j in (5, 8, 11):
        cnt = Counter(phase(a, j) for a in range(1, 1 << j, 2))
        vals = set(cnt.values())
        assert len(vals) == 1, f"j={j}: not equidistributed: {sorted(vals)}"
        print(f"j={j:2d}: phase exactly equidistributed over odd a mod 2^{j} "
              f"({len(cnt)} phases x {vals.pop()} families each)")
    # (2) alternator extremality
    for y in range(2, 12):
        A = (4 ** y - 1) // 3
        assert v2(3 * A + 1) == 2 * y
        for j in range(3, 2 * y + 1):
            assert phase(A, j) == 1, (y, j)
    print("alternators: v2(3*A_y+1) = 2y exact and phase r_j = 1 for all "
          "j <= 2y (y=2..11): deepest fuel at first tick, phase-extremal")
