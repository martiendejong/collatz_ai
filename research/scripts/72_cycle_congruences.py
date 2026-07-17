"""R1466-1490 TRACK B: ANCHOR-NUMERATOR CONGRUENCES (cycle equation census).
A positive cycle with r odd steps, halving word (j_1..j_r), sum j, satisfies
n*D = W, D = 2^j - 3^r, W = sum_i 3^{r-1-i} 2^{J_i} (J_i = partial sums).
D|W is NECESSARY. Meet-in-the-middle exact count of words with W = 0 mod D
for j = ceil(r log2 3) (the only positive-D window with small D).
Question: does W mod D equidistribute (hits ~ C(j-1,r-1)/D, congruences add
NOTHING beyond generic 1/D thinning) or is there an obstruction (hits << expected)?"""
import sys, math
from math import comb
from collections import defaultdict
sys.stdout.reconfigure(encoding="utf-8")
ALPHA = math.log2(3)

def census(r, j):
    D = 2 ** j - 3 ** r
    r1 = r // 2; r2 = r - r1
    # left: words (j_1..j_r1), partial value A = sum_{i<r1} 3^{r1-1-i} 2^{J_i}, total halvings s
    # right: words (j_{r1+1}..j_r), value B likewise with its own local partial sums
    # W = 3^{r2} * A + 2^s * B
    def half(rr, smax):
        out = defaultdict(lambda: defaultdict(int))  # s -> (value mod D) -> count
        # iterate compositions of any total s <= smax into rr parts >= 1
        def rec(parts_left, s, J, val):
            if parts_left == 0:
                out[s][val % D] += 1
                return
            # this part's contribution: 3^{parts_left-1} ... careful: value built as
            # val = sum 3^{rr-1-i} 2^{J_i}; add term for position i = rr - parts_left
            i = rr - parts_left
            term = pow(3, rr - 1 - i, D) * pow(2, J, D)
            for js in range(1, smax - s - (parts_left - 1) + 1):
                rec(parts_left - 1, s + js, J + js, val + term)
        rec(rr, 0, 0, 0)
        return out
    L = half(r1, j - r2)   # left needs >= 1 per remaining right part
    R = half(r2, j - r1)
    hits = 0
    p3 = pow(3, r2, D)
    for s, dl in L.items():
        need_right = j - s
        if need_right < r2: continue
        dr = R.get(need_right)
        if not dr: continue
        p2s = pow(2, s, D)
        # need 3^{r2} A + 2^s B = 0 mod D  ->  B = -3^{r2} A / 2^s mod D
        inv = pow(p2s, -1, D)
        for a, ca in dl.items():
            b_need = (-p3 * a * inv) % D
            cb = dr.get(b_need)
            if cb: hits += ca * cb
    total = comb(j - 1, r - 1)
    return D, total, hits

print(f"{'r':>3s} {'j':>3s} {'D':>14s} {'#words':>14s} {'hits':>6s} {'expected':>10s} {'ratio':>7s}")
for r in range(1, 25):
    j = math.ceil(r * ALPHA)
    if 2 ** j - 3 ** r <= 0: j += 1
    D, total, hits = census(r, j)
    exp = total / D
    ratio = hits / exp if exp > 0 else float('nan')
    star = "  <-- trivial cycle" if r == 1 else ""
    print(f"{r:3d} {j:3d} {D:14d} {total:14d} {hits:6d} {exp:10.3f} {ratio:7.2f}{star}")

print("\nAny hit = word passing the congruence test (necessary, not sufficient:")
print("2-adic valuations must also match). ratio ~ 1 across the board would mean")
print("W mod D equidistributes: congruences give only generic 1/D thinning (no-go).")
print("ratio -> 0 would be a congruence obstruction = a new attack on cycles.")
