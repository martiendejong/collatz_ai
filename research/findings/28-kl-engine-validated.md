# F28 — G1b: the K–L Induction Engine Validated Against the Historical Plateau

**Status: PIPELINE VALIDATED. Certified-form exponents reproduce the pre-2003 state of the art.**

## The two formulations tried
1. **Min-refinement Perron operator: FAILED** (γ collapsed to the bisection floor). Diagnosis: the
   min gives an adversary compounding sabotage power over unknown trits — every branch discounted
   by future worst-case, Perron vector collapses to chains. Documented as a dead end.
2. **K–L induction (correct):** state = (parity, o, r mod 3^(K−o)) — the budget level is always
   K − o, collapsing the state space. Harvest leaf mass Σθ^γ with θ = 3^o/2^e < 1 at
   info-exhaustion (o = K) or horizon D; bisect γ where worst-root mass = 1. Then
   π(x) ≥ c·x^γ follows by induction (π_total ≥ π_worst-class closes the recursion).

## Results
| K (mod 3^K) | D | certified γ ≥ |
|---|---|---|
| 3 | 30 | 0.6012 |
| 4 | 40 | 0.6306 |
| 5 | 50 | 0.6446 |
| 6 | 60/80 | 0.6663 (D-saturated) |
| 7 | 80 | 0.6815 |

Matches the Applegate–Lagarias tree-search era (~0.65–0.68) — exactly where the literature stood
before Krasikov–Lagarias 2003. **The pipeline is historically calibrated.**

## Caveats and next steps (A1 continuation)
- Root sampling at K=7 (200 of ~2916 roots) — full-root sweep needed for theorem-grade output.
- The gap 0.68 → 0.84 is the K–L difference-inequality amortization (system of coupled
  inequalities across classes and scales, LP-solved) — next layer to implement.
- The gap 0.84 → beyond is depth: K = 10–12 with the amortization, modern LP. That is the
  publishable target.

Related: [[27-backward-tree-g1]], PLAN.md Track A1
