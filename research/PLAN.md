# PLAN — The Road From Here (drafted 2026-07-12, after 9 rounds)

Objective: solve the Collatz conjecture, or failing that, produce the maximum amount of *provable*
progress and position the reduction so the community can finish it. Three tracks, decision gates,
honest probabilities.

---

## Track A — Provable increments (compute-assisted theorems). Highest expected value.

### A1 · FLAGSHIP: break the Krasikov–Lagarias 0.84 barrier
The best unconditional density result (≥ x^0.84 of integers reach 1) dates from 2003, was
computer-aided at shallow modulus depth (3-adic systems solved partly by hand), and has not been
improved in 22 years. Modern LP solvers + our exact (a,k) transition algebra can generate and
certify their difference-inequality systems at far deeper moduli (3^7…3^11).
- Step 1: reproduce 0.84 from their published system (pipeline validation).
- Step 2: auto-generate the inequality system at depth 3^k, solve LP, extract exponent.
- Step 3: rational-arithmetic certificate (exact fractions, no floats) → theorem-grade.
- Success = any exponent > 0.84: a publishable unconditional theorem. P(success) ≈ 30–50%.
- Risk: structural ceiling of the method; check follow-up literature first (nothing surfaced 2026).

### A2 · Mechanize m-circuit exclusion past the ~77/91 record
Our chain algebra generates m-circuit Diophantine systems mechanically; add LLL/Baker–Davenport
reduction per shape family. Success = new cycle-exclusion record. P ≈ 20%.

### A3 · The B mod M obstruction (cheap exploration)
In Z/(2^S − 3^K): 2^S ≡ 3^K collapses B(shape) — compute closed forms, test equidistribution of
B mod M over falling shapes, hunt small-prime obstructions. Any structure found feeds A2. P ≈ 10%.

## Track B — The observatory (standing model tests, low cost)
- B1: port the record scanner to C/numba; extend to 10^11–10^12. Every record = a hypothesis test
  of the tail law at new depth. A record OFF the curve would be a major discovery in either direction.
- B2: automate census periods 6–8; publish census tables with the note.

## Track C — The wall itself (theory; where the proof lives)

### C1 · Position and publish the reduction (immediate)
- Literature-diff Theorem C against arXiv 2603.25753 ("one-bit orbit mixing") and Kontorovich–Sinai.
- Polish NOTE.md (now 9 theorems) into an arXiv-format note.
- Submit the elementary theorems to the ccchallenge.org Lean pipeline (all are formalizable).
- Goal: put H(ε)/ε=¼ in front of people with Diophantine-equidistribution tools.

### C2 · Extend Theorem 9 from periodic to structured-aperiodic
Theorem 9 (new, unconditional) kills every periodic drift-defiance mode at logarithmic length via
2-adic separation (n cannot agree with a non-integer 2-adic point beyond log₂n digits).
Program: extend to (a) eventually-periodic patterns (immediate), (b) patterns with periodic
"skeleton" and bounded defects (finite unions — should work), (c) low-complexity aperiodic
patterns (automatic/substitutive streams — their realizations are 2-adic "automatic numbers";
separation from ℕ may be provable via Mahler/Christol theory: **is the realization of an automatic
climbing stream ever a positive integer?** — concrete, novel, attackable).
The gap to close: arbitrary streams have realization sets of full Hausdorff dimension ~0.95 in Z₂ —
entropy is the enemy; each complexity class we kill shrinks the viable conspiracy space.

### C3 · Short-interval unfairness bounds (the endgame shape)
Target theorem: "no orbit at value v is ≥26% unfair for m consecutive macro-steps once
m ≥ f(v)" for ANY explicit f. Then: f ≲ (log v)/(constants) ⟹ full conjecture (via excursion
budget); any f ⟹ quantitative divergence-rate bounds. The needed input is effective
equidistribution of 3^j mod 2^i along orbit sequences — Baker/Ridout/exponential-sum territory.
Theorem 9 is the m = O(log v) statement for periodic windows; C3 asks for it with the period
assumption removed. This is the precise, minimal statement of the whole problem.

## Decision gates
- G1 (next session): A1 step 1 — reproduce 0.84. If pipeline works → A1 becomes the main line.
- G2: C1 literature diff — if Theorem C is subsumed, pivot C1 to contributing Theorem 9 + census
  + tail law to the existing efforts; if novel in form, arXiv the note.
- G3 (after A1 outcome): success → write paper; failure → autopsy feeds A3/C2.
- Standing: B-track runs in background always; every anomaly halts everything else.

## Honesty clause
Tracks A and B produce theorems and instruments with high probability. Track C is the conjecture
itself; no plan makes it likely — the plan makes contact with it *possible* by keeping every
attack aimed at the one named assumption H(ε) and its 2-adic separation structure (Theorem 9's
mechanism is the only unconditional per-orbit tool we have found; C2 grows it).
