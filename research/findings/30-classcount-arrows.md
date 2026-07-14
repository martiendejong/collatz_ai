# F30 — The Fine→Coarse Generation System: sound, tight, closure pending

**Status: inequality family DERIVED + verified sound and 99.95% tight; certificate closure open.**

## Round-14 failures (recorded)
- Multi-level class-count operator with per-level unfolding: **structural dead end** — odd classes
  at the finest modulus have no identity within budget (their rule needs a deeper trit), get forced
  to 0, and poison the cascade through the parity-children. Not a bug; a wrong architecture.

## The correct arrow family (derived, verified)
For odd m ∈ P: 3m+1 = 2^j·n with (j, n) UNIQUE. Conversely every (n ∈ P, valid j) generates
m = (2^j·n − 1)/3 ∈ P — unique source per target, hence no collisions and:

> π_{c mod 3^(k−1)}(x) ≥ Σ_{(a mod 3^k, j) → c} π_{a}^{odd}(3x/2^j)

Arrows flow FINE → COARSE (one trit consumed per generation — the 2-adic-in/3-adic-out entropy
core made visible), scales 3/2^j < 1 for j ≥ 2. Numerical soundness check at k=3, x=2×10⁵:
lhs ≥ rhs in all 9 target classes, with rhs = 99.95% of lhs — the family is essentially an
EQUALITY (as the uniqueness argument predicts). This is the legitimate class-count backbone.

## The remaining puzzle: closure
The arrows only push mass downward in modulus; the finest level needs feeding. Pigeonhole
(some refinement ≥ coarse/3) fails per-class (which refinement is unknown). This closure is
exactly what the Krasikov–Lagarias LP resolves — their paper's actual system is the missing piece.

## Next round (A1)
Fetch Krasikov–Lagarias 2003 (arXiv math/0205002, "Bounds for the 3x+1 problem using difference
inequalities"), implement their exact system on our verified arrow family, reproduce 0.84,
then push depth.

Related: [[29-chainsafe-optstop]], PLAN.md Track A1
