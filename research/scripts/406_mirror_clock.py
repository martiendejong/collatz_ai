"""Obs 615: LEMMA 38b -- the mirror clock dichotomy (proof + verification).

Statement. For j >= 3 the unit group mod 2^j splits as <3> u (-<3>)
(index 2, and -1 is not in <3> since -1 = 7 mod 8 while <3> = {1,3}
mod 8). Hence for every odd a exactly one of the two clocks resonates:
  a = 1,3 mod 8  (a in  <3>-coset): MIN clock v2(a*3^k - 1) is an
      affine tower (Lemma 38); PLUS clock capped at v2 <= 2.
  a = 5,7 mod 8  (a in -<3>-coset): PLUS clock v2(a*3^k + 1) is an
      affine tower; MIN clock capped (Lemma 38).
Proof of the caps: depth >= 3 on the min clock needs 3^k = a^(-1)
mod 8 with a^(-1) in {5,7} -- impossible; on the plus clock it needs
3^k = -a^(-1) mod 8 with -a^(-1) in {5,7} -- impossible. Towers by
cyclicity of <3> mod 2^j, exactly as in Lemma 38.

Consequences.
  C1: deep v2(a*3^k + 1) = deep trailing ONES of a*3^k = fuel. The
      rich refill caste {5,7 mod 8} is exactly the plus-resonant coset:
      caste = coset character. That character is the Kronecker symbol
      (-2/a): +1 iff a = 1,3 mod 8.
  C2: alternators (4^y-1)/3 = 5 mod 8 (y >= 2), i.e. maximally
      plus-resonant -- the gateways of 1 are fuel-clock families.
  C3: unified clock W_a(k) = max of the two: EVERY family carries
      exactly one affine tower; the sign is the caste character.
"""

def v2(x: int) -> int:
    v = 0
    while x % 2 == 0:
        x //= 2
        v += 1
    return v

def tower_check(vals, K):
    """{k : val >= j} must be one residue class mod 2^(j-2), nested."""
    phases = []
    for j in range(3, 10):
        hits = [k for k, x in enumerate(vals, start=1) if x >= j]
        mod = 2 ** (j - 2)
        assert hits, f"no depth-{j} hits"
        cls = {k % mod for k in hits}
        assert len(cls) == 1, f"j={j}: {len(cls)} classes"
        assert abs(len(hits) - K // mod) <= 1, f"j={j}: sparse"
        phases.append(cls.pop())
    for (j, r), rn in zip(enumerate(phases, start=3), phases[1:]):
        assert rn % (2 ** (j - 2)) == r, f"not nested at j={j}"
    return phases

if __name__ == "__main__":
    K = 1024
    n_min = n_plus = 0
    for a in range(1, 402, 2):
        minus = [v2(a * 3 ** k - 1) for k in range(1, K + 1)]
        plus = [v2(a * 3 ** k + 1) for k in range(1, K + 1)]
        if a % 8 in (1, 3):
            assert max(plus) <= 2, f"a={a}: plus clock not capped"
            tower_check(minus, K)
            n_min += 1
        else:
            assert max(minus) <= 2, f"a={a}: min clock not capped"
            tower_check(plus, K)
            n_plus += 1
    print(f"Lemma 38b verified: {n_min} min-resonant (a=1,3 mod 8), "
          f"{n_plus} plus-resonant (a=5,7 mod 8), all towers affine+nested")
    # C2: alternators
    for y in range(2, 20):
        assert ((4 ** y - 1) // 3) % 8 == 5
    print("C2 verified: (4^y-1)/3 = 5 mod 8 for y=2..19 (plus-resonant coset)")
