"""Obs 611: Rule 30 lens, tasks S1+S2.

S1 (flip-parity criterion): multiply-by-3 on binary digits is a carry
automaton; bit j of 3^k evolves as
    bit_j(3^(k+1)) = bit_j(3^k) XOR driver_j(k),
    driver_j(k) = bit_{j-1}(3^k) XOR carry_j(k),
where the driver is a function of the lower block 3^k mod 2^j only,
hence periodic in k with the lower-block period p_j = ord(3 mod 2^j)
= 2^(j-2) (j >= 3). Rule 30's doubling criterion: the period of bit j
doubles (relative to the lower block) iff the number of driver hits per
lower period is ODD; if EVEN, bit j is a branch point (seed information
survives). Claim: for the Collatz clock the parity is odd at EVERY
level (equivalent to LTE ord(3 mod 2^(j+1)) = 2^(j-1)) -- the ordered
side of Collatz has no branch points at all.

S2 (family clocks): for odd a, the sequence B_a(k) = v2(a*3^k - 1).
Prediction from subgroup structure (<3> mod 2^j = residues {1,3} mod 8,
index 2): a == 1,3 (mod 8)  -> full ruler, unbounded, +1 per doubling
                               window, deep set {v2 >= j} = one residue
                               class mod 2^(j-2) (a phase);
            a == 5,7 (mod 8) -> capped clock, v2 <= 2 forever (flat).
Then all family-dependence sits in the PHASE WORD (which lift succeeds
at each level) = the binary digits of the 2-adic discrete log
log_3(a^{-1}), and the innovation VALUES are affine for every family.
"""

def v2(x: int) -> int:
    v = 0
    while x % 2 == 0:
        x //= 2
        v += 1
    return v

# ---------- S1: flip parity at every clock level ----------
def flip_parity(j: int) -> int:
    """Parity of the number of k in one lower period where bit j flips."""
    p = 2 ** (j - 2)  # period of 3^k mod 2^j, j >= 3
    m = 2 ** (j + 1)
    x = 1
    flips = 0
    for _ in range(p):
        y = (3 * x) % m
        if ((x >> j) ^ (y >> j)) & 1:
            flips += 1
        x = y
    return flips & 1

def check_S1(jmax: int = 20):
    pars = [flip_parity(j) for j in range(3, jmax + 1)]
    assert all(p == 1 for p in pars), f"even flip count found: {pars}"
    return pars

# ---------- S2: family clocks ----------
def family_seq(a: int, K: int):
    return [v2(a * 3 ** k - 1) for k in range(1, K + 1)]

def check_S2(K: int = 1024, amax: int = 201):
    flat, ruler = [], []
    for a in range(1, amax + 1, 2):
        s = family_seq(a, K)
        if a % 8 in (5, 7):
            assert max(s) <= 2, f"a={a}: capped-clock prediction fails, max={max(s)}"
            flat.append(a)
        else:
            # affine tower: {k : v2 >= j} is exactly ONE residue class mod
            # 2^(j-2), full within range, nested with index exactly 2.
            # (Window-max is NOT monotone +1: phases can place deep hits
            # early, e.g. a=11 has v2=5 at k=1 -- that is the phase effect,
            # not a violation.)
            phases = []
            for j in range(3, 10):
                hits = [k for k, x in enumerate(s, start=1) if x >= j]
                mod = 2 ** (j - 2)
                assert hits, f"a={a}: no depth-{j} hits in range"
                cls = {k % mod for k in hits}
                assert len(cls) == 1, f"a={a}, j={j}: {len(cls)} classes"
                r = cls.pop()
                assert abs(len(hits) - K // mod) <= 1, f"a={a}, j={j}: sparse class"
                phases.append(r)
            for (j, r), r_next in zip(enumerate(phases, start=3), phases[1:]):
                assert r_next % (2 ** (j - 2)) == r, f"a={a}: tower not nested at j={j}"
            phase_word = [
                (r_next - r) // (2 ** (j - 2))
                for (j, r), r_next in zip(enumerate(phases, start=3), phases[1:])
            ]
            ruler.append((a, phase_word))
    return flat, ruler

if __name__ == "__main__":
    pars = check_S1(20)
    print(f"S1: flip parity ODD at every level j=3..20: {pars} -> no branch points")
    flat, ruler = check_S2()
    print(f"S2: a=5,7 mod 8: {len(flat)} families ALL capped at v2<=2 (flat clock)")
    print(f"    a=1,3 mod 8: {len(ruler)} families ALL affine towers: one class per")
    print("    depth, nested index 2, never capped. Example phase words (j=3->9):")
    for a, pw in ruler[:8]:
        print(f"      a={a:3d}: {pw}")
    print("OK: innovation values affine for every family; family identity lives")
    print("    only in the phase word (2-adic discrete log).")
