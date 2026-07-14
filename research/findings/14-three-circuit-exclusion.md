# F14 — Three-Circuit Cycles: none with any segment index ≤ 24

**Status: PROVEN (exhaustive in range: 36,019,816 shapes, exact integer arithmetic).**

Chain composition of the exact transition gives, for a 3-circuit (a₁,k₁)→(a₂,k₂)→(a₃,k₃)→(a₁,k₁):
a₁·(2^(K+W) − 3^K) = B(shape), with B accumulated iteratively. Pruning: 2^(K+W) must sit within
8 bits above 3^K (else a₁ < 1 or non-integral).

**Result: only the trivial cycle (all segments (1,1), a=1).** Elements were UNBOUNDED —
only segment geometry (kᵢ ≤ 24, band 8) was constrained. Combined ledger: 1-circuits none (Steiner,
range k,w<400), 2-circuits none (shapes ≤ 60), 3-circuits none (kᵢ ≤ 24), m ≤ ~77 none (literature),
length ≥ 1.69e11 (F8).

Related: [[13-two-circuit-exclusion]], [[08-cycles-diophantine]]
