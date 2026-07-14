# F27 — G1/G2: Backward-Tree Engine v0, the Correct Functional, and the Theorem C Diff

**Status: G1 partial (machinery validated, functional identified); G2 resolved.**

## G2: Theorem C vs Chang (arXiv 2603.25753)
Chang reduces Collatz to EXACT one-bit balance (bit 4 at burst-ending times, two classes mod 32,
sparse subsequence). Structurally different from Theorem C: his hypothesis is sharp and fragile
(exact balance), ours is crude and robust (25% tolerance on graded frequencies, both halves at once,
elementary proof). **Not subsumed — complementary. NOTE.md remains worth publishing.**

## G1: backward-tree engine v0
State machine (parity, ternary residue, budget) with transitions
m → 2m (free) and m → (m−1)/3 (iff m even, ≡1 mod 3; consumes one trit) — validated.
Worst-case uniform-depth min-leaf DP up to (K=8, D=32): γ ≥ 0.3237, growing slowly in K, D.

## The lesson that matters (the round's real yield)
Uniform-depth counting is the WRONG functional: backward /3-steps DECREASE values, so depth-D
leaves are not comparable in size — small leaves pack extra predecessors below x that depth
counting ignores. The correct object, derived from our own recursion:

> |P(x)| = |P(x/2)| + |P^(≡2 mod 3)(1.5x)| — a TWO-SCALE identity, refined over classes mod 3^K.

The exponential ansatz π_r(x) = u_r·x^γ turns this into the nonlinear system
u_r = 2^(−γ)·u_{c1(r)} + 1.5^(γ)·u_{c2(r)} — precisely the Krasikov–Lagarias difference-inequality
structure (their rigor comes from careful truncation of the upward scale 1.5x). Pipeline for next
session: (1) build the class-coupled system mod 3^K, (2) solve the eigenproblem for heuristic γ(K),
(3) implement the K–L truncation to make it rigorous, (4) push K beyond their 2003 depth.

Related: [[24-theorem-c-conditional]], [[19-realization-census]], PLAN.md Track A1
