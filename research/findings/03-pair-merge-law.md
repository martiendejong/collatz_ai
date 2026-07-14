# F3 — The Pair-Merge Law: exact, exhaustively verified, and unique to multiplier 3

**Status: PROVEN mechanism + census 3,002/3,002 (100%).**

## Mechanism
Consecutive sequences (k, k+1) in family a have segment endpoints x = a·3^k − 1 and x′ = a·3^(k+1) − 1 satisfying **x′ = 3x + 2**. Halving both:

> x′/2 = 3·(x/2) + 1 = exactly the Collatz step applied to y = x/2, whenever y is odd.

So the two members of a "true pair" (those with y odd) merge deterministically just after their segments end.

## Census (script 03)
- **3,002 of 3,002** true pairs across families a < 2000 merged. Zero exceptions.
- Family pairing offset (whether pairs are (1,2),(3,4)… or (2,3),(4,5)…) splits **333 / 334** — exactly half-half, decided by the parity of y in each family.

## Uniqueness to c = 3 (the near-miss theorem)
For a general map cn+1 the identity becomes x′/2 = c·y + (c−1)/2 while the map needs c·y + 1. Equal **iff (c−1)/2 = 1 iff c = 3.**
For 5n+1 the mismatch is a permanent off-by-one; measured would-be merge points: 187 vs 186, 937 vs 936, 2187 vs 2186. Only 18/199 of 5n+1 (n,4n+3) pairs merge (coincidence level) vs 100% law-driven in Collatz.

## The 3n−1 warning (what the law cannot prove)
3n−1 has the **same** pair-merge law (x′/2 = 3y − 1 = its own map step): 0/600 true pairs split across cycles. Yet 3n−1 has **three cycles** (min elements 1, 5, 17) and starts 1–10,000 split 3,244/3,213/3,543 across them. **The pair law merges within a tree; it cannot count trees.**

## Litmus doctrine (adopted for every candidate argument)
1. Must fail for 5n+1 (else it ignores the drift).
2. Must fail for 3n−1 (else it disproves cycles that actually exist).
3. Must not abolish the trivial cycle {1,4,2}.

Related: [[04-funnel-and-controls]], [[08-cycles-diophantine]]
