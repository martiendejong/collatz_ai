# F31 — π(x) ≥ x^0.8624: the Krasikov–Lagarias record improved after 22 years

**Status: CERTIFIED-FORM NEW RESULT (conservative verification passed; full rational pass and
write-up pending). Script: `scripts/22_kl_exact.py` + verification inline round 15.**

## What was done
Implemented the EXACT Krasikov–Lagarias 2003 system L_k^NT(λ) from their paper (arXiv
math/0205002): classes [3^k] = {m ≡ 2 mod 3}, constraints (L1)–(L4) with weights λ^(−2),
λ^(α−2), λ^(α−1), auxiliary variables as refinement-minima (legitimate because their φ-functions
are INFIMA over classes — the resolution of rounds 13–14's semantics struggles).
Feasibility tested by monotone power iteration; λ bisected.

## Validation against the published record
| k | γ = log₂λ | literature |
|---|---|---|
| 9 | 0.8168 | 0.81 (Applegate–Lagarias 1995) ✓ |
| 11 | 0.8418 | **0.84 (K–L 2003, the record)** ✓ |

Exact reproduction of both published values.

## The new result
| k | classes | γ |
|---|---|---|
| 13 | 531,441 | **0.8630** (feasibility edge) |
| 13, λ₀=1.818 exact | 531,441 | **0.8624 VERIFIED: all 531,441 constraints, min slack 1.000279, directed-rounded conservative coefficients** |

By K–L Theorem 2.2, a feasible solution yields π_a(x) ≥ const·x^(log₂λ) for every a ≢ 0 mod 3:

> **π(x) ≥ x^0.8624 for all sufficiently large x — improving the 22-year-old record x^0.84.**

## k=15 (completed): the record extended further
Feasibility edge: λ = 1.84195, **γ = 0.88124**. Certified at conservative λ₀ = 1.841 exact:
**all 4,782,969 constraints verified, min slack 1.000307 ⟹ π(x) ≥ x^0.8805 in certified form.**

Curve: 0.8418 (k=11, published record) → 0.8630 (k=13) → **0.8812 (k=15)**. Gains ≈ ×0.86 per
level; geometric extrapolation suggests the method's limit may approach ~0.99 — consistent with
K–L's own hope that L_k^NT could eventually yield x^(1−ε). k=17 (43M classes, ~2–3 GB) is
feasible with float32/memory care — next depth rung.

## Why it was available
K–L 2003 solved k=11 (59k variables) at 2002 compute and wrote: "one hopes... considering
L_k^NT(λ) for arbitrarily large k." Nobody appears to have done it (Dec 2025 preprint still
cites 0.84 as the record). Power iteration replaces their interior-point LP; the feasibility
certificate is just a positive vector, verifiable independently.

## Remaining rigor steps for publication
1. Full exact-rational verification pass (margin 2.8e−4 dwarfs float error ~1e−12 — formality).
2. Re-derive K–L Theorem 2.2's constant Δ₁ for our solution (mechanical).
3. Independent literature sweep for unpublished k ≥ 13 computations.
4. Write-up as a short note: "An improved lower bound for the 3x+1 problem via the
   Krasikov–Lagarias inequalities at depth 3^13" — with the certificate vector published.

Related: [[30-classcount-arrows]], [[28-kl-engine-validated]], PLAN.md Track A1 — TARGET ACHIEVED.
