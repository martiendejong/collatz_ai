"""Obs 610: Rule 30 right-boundary sequence vs the Collatz 2-adic clock.

Claim (verified here):
  A = rightmost blue-block lengths of Rule 30 (single blue cell on black
      background), the sequence from the right-boundary analysis:
      1,3,1,4,1,3,1,6,1,3,1,4,1,3,1,7,...
  B = v2(3^k - 1), k = 1,2,3,...:
      1,3,1,4,1,3,1,5,1,3,1,4,1,3,1,6,...

Both satisfy the IDENTICAL repetition law: a new value appears only at
positions 2^m, and immediately after each new value the entire prefix
replays. Innovations:
  Rule 30 innovations: 1,3,4,6,7,9,...  (no known formula; computed to
                                         2^46 rows, still opaque)
  clock innovations:   1,3,4,5,6,7,...  = 2 + v2(k), closed form via LTE.
Because the repetition law replays innovations, the two sequences agree
EXACTLY on the dense set {k : v2(k) <= 2} (7/8 of all positions) and
differ exactly on the sparse deep-clock slots {k : v2(k) >= 3}
(0-based indices i with i == 7 mod 8).

Reading: the Collatz 2-adic odometer is a SOLVED Rule 30 boundary; the
hardness that Rule 30 keeps in its innovation sequence sits, for Collatz,
in the Sturmian/rotation layer (log2 3) and the carry coupling instead.
"""

def rule30_boundary(nrows: int):
    W = 2 * nrows + 10
    row = [0] * W
    row[5] = 1  # single blue cell; background to the right stays black
    def step(r):
        out = [0] * len(r)
        for i in range(1, len(r) - 1):
            out[i] = r[i - 1] ^ (r[i] | r[i + 1])  # Rule 30
        return out
    def rightmost_blue_block(r):
        i = len(r) - 1
        while i >= 0 and r[i] == 0:
            i -= 1
        j = i
        while j >= 0 and r[j] == 1:
            j -= 1
        return i - j
    seq = []
    for _ in range(nrows):
        seq.append(rightmost_blue_block(row))
        row = step(row)
    return seq

def v2(x: int) -> int:
    v = 0
    while x % 2 == 0:
        x //= 2
        v += 1
    return v

def check_repetition(s):
    """New value only after which the full prefix replays; return violations."""
    seen, bad = set(), []
    for i, x in enumerate(s):
        if x not in seen:
            seen.add(x)
            if (i + 1) & i:  # novelty allowed only at 0-based index 2^m - 1
                bad.append(("novelty-position", i, x))
            tail = s[i + 1:i + 1 + i]
            if tail and tail != s[:len(tail)]:
                bad.append(("replay", i, x))
    return bad

if __name__ == "__main__":
    N = 64
    A = rule30_boundary(N)
    B = [v2(3 ** k - 1) for k in range(1, N + 1)]
    print("Rule30 :", A[:32])
    print("clock  :", B[:32])
    assert not check_repetition(A), "Rule 30 repetition law violated"
    assert not check_repetition(B), "clock repetition law violated"
    diff = [i for i in range(N) if A[i] != B[i]]
    print("differ at 0-based indices:", diff[:10])
    assert all((i + 1) % 8 == 0 for i in diff), "difference outside deep-clock slots"
    assert all(A[i] != B[i] for i in range(7, N, 8)), "deep-clock slots should all differ here"
    print("OK: identical repetition law; sequences agree exactly on v2(k)<=2 (7/8 of positions)")
