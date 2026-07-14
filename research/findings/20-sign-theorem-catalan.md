# F20 — The Sign Theorem and the Catalan Origin of the Known Cycles

**Status: PROVEN.**

## The Sign Theorem
For any cyclic shape, the closure constant B(shape) = Σᵢ (2^wᵢ − 1)·3^(k's after i)·2^(S's before i)
is a sum of positive terms, hence **B > 0 always** (verified 20,000 random shapes; obvious from the form).
Therefore, with a₁ = B/(2^(K+W) − 3^K):

> **sign(a₁) = sign(2^(K+W) − 3^K).**
> A POSITIVE cycle requires a net-falling shape (W > K·(log₂3 − 1));
> every net-climbing shape realizes on the NEGATIVE side.

This is why all the "conspiracy" objects (the 9/8 ladder limit −5, the −17 cycle) are negative:
climbing shapes are structurally banished from the positive integers. The divergence question is
what remains: climbing *aperiodic* streams (no closure constraint) — E★ territory.

## The Catalan origin of the known cycles
Integer realization needs (2^S − 3^K) | B. The gap table (min over S for each K):
K=1: **+1** (2²−3) · K=2: **−1** (2³−9) · K=3: +5 · K=4: −17 · K=7: −139 · … growing forever.

|gap| = 1 makes integrality FREE — and by **Mihailescu's theorem (Catalan's conjecture)**, 8 and 9
are the only consecutive proper powers: the unit gaps at (K,S) = (1,2) and (2,3) are the only ones
that will ever exist. Hence:

> **The trivial cycle (n=1) and the −5 cycle are exactly the two Catalan solutions.**
> n=1: gap +1 (net-falling, positive side). n=−5: gap −1 (net-climbing, negative side).
> Every further cycle (like −17, gap −139) needs a divisibility coincidence B ≡ 0 mod |gap|,
> and Baker-type bounds make |gap| grow superpolynomially — the throttle behind all length bounds.

## Consequences for the conjecture
The cycle half now reads: no net-falling shape beyond (1,1) has (2^S − 3^K) | B(shape).
The two structural "gifts" (unit gaps) are both spent — one on the trivial cycle, one on −5.
Arithmetic has no more free cycles to give; everything else must be a coincidence, and the
coincidences are Baker-throttled. This is, in compressed form, WHY the conjecture is plausible:
**the only cheap cycles that could exist, do exist, and they are 1 and −5.**

Related: [[19-realization-census]], [[08-cycles-diophantine]], [[17-ladders-and-minus5]]
