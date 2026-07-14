# F13 — Two-Circuit Cycles: none with any segment shape up to 60

**Status: PROVEN (exhaustive in range).**

## Setup
A 2-circuit cycle is a chain (a₁,k₁) → (a₂,k₂) → (a₁,k₁) of two macro-steps. Composing the exact
transition algebra (F1) and clearing denominators gives the linear condition

> a₁ · (2^(k₁+k₂+w₁+w₂) − 3^(k₁+k₂)) = 3^(k₂)·(2^(w₁) − 1) + (2^(w₂) − 1)·2^(w₁+k₂)

## Search
All (k₁, w₁, k₂, w₂) with every component ≤ 60 — 12.9 million shapes, exact integer arithmetic.
**Result: the only solution is (1,1,1,1) = the trivial cycle traversed twice.**

## Context
Consistent with and complementary to Simons–de Weger (no m-circuit cycles for m ≤ ~77, by
transcendence bounds). Our range covers a different slice: arbitrary a₁ (unbounded!) with segment
shapes ≤ 60 — the transcendence results bound m, this bounds the per-segment geometry with NO
bound on the elements. The two exclusion regions overlap but neither contains the other.

## Open frontier refined
A nontrivial cycle now requires: ≥ 1.69e11 total steps (F8), > ~77 circuits (literature),
AND at least one segment parameter > 60 (this result). The (a,k) chain formulation makes
m-circuit systems mechanical to generate; pushing B and m further is pure compute.

Related: [[08-cycles-diophantine]], [[01-family-pair-coordinates]]
