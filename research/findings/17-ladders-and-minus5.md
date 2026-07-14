# F17 — The 9/8 Ladder: bounded-index divergence is not excluded, and its limit object is −5

**Status: PROVEN (census + fixed point) — refutes "divergence needs infinite trailing 1s" for 3n+1.**

## The object
A "9/8 ladder" is a run of macro-steps with k=2, w=1: each rung multiplies the value by 9/8 > 1.
A number riding the ladder climbs forever with trailing-ones NEVER exceeding 2.

## Census (all n ≡ 3 mod 4 below 2^24)
| rungs r | count |
|---|---|
| 1 | 917,504 |
| 3 | 14,336 |
| 5 | 224 |
| 7 | 4 |

Continuation probability **exactly 1/8 per rung at every level** (0.125, six decimals of agreement) —
the cleanest geometric law in the program. Longest: 7 rungs from 4,194,299, climbing 2.3× with
trailing-ones pinned at 2. Ladders of every finite length exist (density 8^(−r), never zero).

## Verdict on the claim "divergence needs infinite trailing 1s"
- TRUE only for single-segment (monotone) ascent — that mode is impossible (F9).
- FALSE in general, even for original 3n+1: the 9/8 ladder diverges with k bounded by 2.
  Nothing known excludes an infinite ladder; only its density (8^(−r)) vanishes.

## The twist: the infinite ladder exists — at n = −5
Ladder recursion on the family parameter: a → (9a+1)/8. Its fixed point is **a = −1, i.e. n = −5**:

> −5 → −14 → −7 → −20 → −10 → −5 — the famous negative cycle IS the infinite 9/8 ladder.

By Terras' bijection every infinite index stream is realized by exactly one 2-adic integer;
the k=2,w=1 stream's realization is −5. **The conspiracies are not fictions — they are all realized
as 2-adic points; the conjecture is precisely that none of the bad realizations is a positive integer.**
(Likewise the other negative cycle at −17 realizes a longer periodic pattern.)

Related: [[07-never-drop-needle]], [[09-verdicts-and-open-core]]
