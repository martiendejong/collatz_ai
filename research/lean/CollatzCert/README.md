# CollatzCert — a Lean-verified feasibility certificate beating the Krasikov–Lagarias record

This is a **mathlib-free Lean 4 project** that machine-verifies, in exact
integer arithmetic, a feasibility certificate for the Krasikov–Lagarias
difference-inequality system `L_12^NT(λ)` at

```
γ = 213/250 = 0.852        λ = 2^γ
```

By **Theorem 2.2 of Krasikov & Lagarias**, *Bounds for the 3x+1 Problem using
Difference Inequalities*, Acta Arith. 109 (2003) 237–258 (peer-reviewed,
not formalized here), feasibility of `L_k^NT(λ)` implies

```
π₁(x) > x^0.852   for all sufficiently large x,
```

where π₁(x) counts integers below x whose 3x+1 orbit reaches 1. The published
record exponent is **x^0.84** (k = 11, loc. cit.; confirmed as the standing
record by arXiv:2512.13760, Dec 2025). This certificate therefore witnesses a
strict improvement — with the witness checking done by Lean, not by trust.

## What Lean verifies (theorems in `CollatzCert.lean`)

1. `coefA`, `coefB1`, `coefB3` — the dyadic coefficients `pa/2^Q ≤ λ⁻²`,
   `p1/2^Q ≤ λ^(log₂3−2)`, `p3/2^Q ≤ λ^(log₂3−1)` via **pure integer power
   inequalities** (e.g. `pa^250 ≤ 2^(250·Q−426)`), which is exactly the
   statement that the dyadics lower-bound the real algebraic coefficients.
2. `vector_ok` — the 177,147-entry certificate vector is positive.
3. `certificate_feasible` — all 177,147 per-class inequalities of
   `L_12^NT(λ)` hold for the vector with the lowered coefficients
   (monotonicity in the coefficients then gives feasibility at the true
   values — a one-line remark, see the docstring).

All proofs are by `native_decide` (compiled kernel-checked evaluation).

## What is cited, not formalized

The reduction *feasibility ⟹ density bound* is Theorem 2.2 of the published
Acta Arithmetica paper. Formalizing that paper is a separate project
(a natural fit for [ccchallenge.org](https://ccchallenge.org/)).

## Build

```
elan default leanprover/lean4:v4.15.0   # or rely on lean-toolchain
lake build
```

## Provenance

Certificate vector produced by power iteration + directed-rounding
verification (scripts 163–166 in `../../scripts/`), exported and
exact-pre-verified in Python big-int arithmetic by `183_lean_export.py`
(margin: min RHS/LHS = 1.000497 over all 177,147 inequalities).
Higher levels are certified numerically up to k = 18 (γ = 0.902); Lean
export of larger k is mechanical (larger data files).

Part of the public research log at
[github.com/martiendejong/collatz_ai](https://github.com/martiendejong/collatz_ai).
