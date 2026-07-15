# OPEN PROBLEMS -- Collatz Research Corpus (compiled 2026-07-15)

Every unsolved thread mentioned anywhere in NOTE.md, REPORT.html (Rounds 1-594),
findings/01-41, PLAN.md, NOTE_DENSITY.md, VERIFICATION.md and COMMUNITY.md.
Ranked: COMPUTABLE first, then LOCAL-ANALYTIC, then HARD.

Difficulty legend:
- COMPUTABLE      -- a finite computation (possibly large) would settle it
- LOCAL-ANALYTIC  -- provable with existing tools, weeks of work
- HARD            -- needs new mathematics

Status snapshot for context: density record pi(x) >= x^0.9146 (k=20, exact-verified,
R576); cycle census complete through period 12 (Thm 17); circuits excluded m <= 6
in-house, m <= 91 by Hercher; tail law tested to 5e10; the conjecture itself OPEN.

---

## COMPUTABLE

### 1. The k=21 Krasikov-Lagarias computation (the standing referee)
- Statement: run L_21^NT (3^21 = 10.5e9 classes, ~16 GB, cloud-scale), measure
  gamma, alpha_21, theta, (a,c); predictions on record: gamma ~ 0.918,
  (a,c) ~ (0.465, 0.528), theta ~ 0.850.
- Where: NOTE.md Thm 16 (falsifiable k=21 predictions); findings/34 (ceiling
  tension, "decisive test: k=21"); findings/40 (R374-380 prediction); REPORT
  R79-85, R374-380, final synthesis ("with the k=21 computation as referee").
- Class: COMPUTABLE.
- Unlocks: strongest out-of-sample test of the lattice model + tempering law;
  new density record x^0.918; primary evidence for item 2.

### 2. The ceiling tension: does q -> 1 or q -> ~0.993?
- Statement: two extrapolations of the min-loss series disagree -- Shanks gives
  q_inf ~ 0.9927 (method ceiling gamma_inf ~ 0.976-0.987), the CV-route gives
  q -> 1 (no ceiling below gamma = 1); flagged "OPEN, DECIDABLE".
- Where: findings/34 (Campaign II, "The ceiling tension"); findings/33 (R39
  gamma-ceiling 0.985-0.987); REPORT R79-85; NOTE.md Thm 12 (q_inf ~ 0.993).
- Class: COMPUTABLE (k=21, then k=23 if still ambiguous, plus longer q-series).
- Unlocks: whether the K-L method can in principle reach x^(1-eps), i.e. whether
  the whole density endgame (items 18-23) is aimed at an attainable target.

### 3. Publication closure of the density record
- Statement: finish the note -- Delta_1 constants for k=17/19/20 (mechanical,
  done for 13/15), refresh stale NOTE_DENSITY.md/tex (still titled "through
  3^19", still lists k=19 as staged; record is now 0.9146 at k=20), third
  literature sweep, arXiv submission, contact Lagarias (Level 6 of the
  verification protocol).
- Where: NOTE_DENSITY.md status checkboxes ([ ] k=19, [ ] contact Lagarias);
  findings/31 ("remaining rigor steps for publication" 1-4); VERIFICATION.md
  Level 6; findings/33 (R31/33 Delta_1 for 13/15 only).
- Class: COMPUTABLE (mechanical + administrative).
- Unlocks: the record becomes a citable theorem; community verification.

### 4. Profile-integrated source law (the queued refinement)
- Statement: the quasi-static law lambda_k = C * CV_top(k) holds only to first
  order (C drifts 5.36 -> 6.36); compute the exact profile-integrated source.
- Where: REPORT R567-575 ("the refinement queued").
- Class: COMPUTABLE (measurement/fit on existing certificates).
- Unlocks: the exact source term needed by the injection-decay statement
  (item 20) and the (a,c) derivation (item 21).

### 5. Exact identification of the correction field g
- Statement: g ~ -log(roulette) explains R^2 = 0.94; identify the residual
  exactly -- spring-blocked first rungs (r = +0.52) and the low-conductor
  character concentration (order-3 character = 28% of AC power, top-3 = 43%).
- Where: findings/39 (R334-340, "the sharpest open experimental question");
  findings/40 (R381-408); REPORT R510-515 (character space).
- Class: COMPUTABLE (finite character/regression analysis on cert_k13..k20).
- Unlocks: a closed form for g turns the spectral-contraction question
  (items 20-23) into a computation about explicit low-order character modes.

### 6. Cycle census extension: period >= 13 and wider symbol bounds
- Statement: extend the exhaustive integer-cycle census past period 12 (current
  bounds: p12 all compositions; p1 k,w <= 40; p2 <= 18; p3 <= 10; p4 <= 7;
  circuits m=5 k <= 6, m=6 k <= 5, m=7 shapes 74.4M).
- Where: NOTE.md Thm 17; findings/19/21/22/32/34; REPORT R421-432; PLAN B2.
- Class: COMPUTABLE (the census program is "mechanical for any bound").
- Unlocks: cycle-freeness complete in ever larger finite regions; complements
  Hercher's m <= 91; publishable census tables.

### 7. Mechanized m-circuit exclusion past Hercher's m <= 91
- Statement: use the exact (a,k) chain algebra to generate m-circuit Diophantine
  systems and kill them per shape family with LLL / Baker-Davenport reduction.
- Where: PLAN Track A2 (P ~ 20%); findings/8 ("open frontier"); findings/13
  ("pushing B and m further is pure compute"); COMMUNITY.md (Hercher 2023).
- Class: COMPUTABLE (standard reduction tooling per family; large but finite).
- Unlocks: a new published cycle-exclusion record; shrinks the only region
  (K > 6.5e10, m > 91) where a nontrivial cycle could live.

### 8. Records / tail-law scan to 1e11-1e12
- Statement: port the sigma-record scanner to C/numba and extend past the
  current 5e10; every new record is a hypothesis test of the exact tail law
  P(sigma >= t) = C * t^(-3/2) * 2^(-t/20) at unprecedented depth.
- Where: PLAN Track B1; findings/26 (status note: "the deepest available probes
  of E-star"); findings/33 (R40 complete to 8e9); REPORT (tested to 5e10).
- Class: COMPUTABLE.
- Unlocks: an on-curve record extends the model's validated range; an OFF-curve
  record "would be a major discovery in either direction".

### 9. Backward-wavefront frontier past 2^25
- Statement: compute exact d(n) beyond 2^25 (next expected seed 63728127 at
  2^26), track the laggard-edge constant c* (17.75 and rising toward the
  stochastic constant ~29).
- Where: findings/39 (R341-353); findings/41 (R433-444, laggards to 2^25).
- Class: COMPUTABLE.
- Unlocks: data on the exact reformulation "bounded c* pointwise <=> conjecture"
  (item 36); locates extremal candidates (ladder-riders in small families).

### 10. Termination-prover search on the Law-A anchor/grammar encoding
- Statement: run matrix/arctic-interpretation (SAT) termination searches on OUR
  anchor/grammar encoding of the rewriting system instead of the generic
  mixed-base encoding used by Yolcu-Aaronson et al. -- "the unopened door".
- Where: REPORT ~R102-103 ("concrete open move"); NOTE.md Thm 14 (gated /
  relative formulations are the legal targets).
- Class: COMPUTABLE (finite SAT search; failure is informative, success huge).
- Unlocks: a machine-found relative-termination certificate for the borrow
  chain -- a genuinely new kind of partial result.

### 11. Certified density exponents for 3n-1 and 5n+1
- Statement: port the K-L engine to the sibling maps and certify the first-ever
  density exponents there (measured: 3n-1 trees ~0.90 each; 5n+1 grounded
  skeletons ~0.55).
- Where: REPORT next-research menu (R119ff); findings/35 (sibling exponents,
  first measurements anywhere).
- Class: COMPUTABLE.
- Unlocks: rigorous grounded-vs-divergent exponent classes; the strongest
  control-system calibration of the whole density method.

### 12. Boundedness test of the alpha = 1.2 max-excursion (open lead #2)
- Statement: test on much wider samples whether the max height-increase at
  alpha ~ 1.2 is bounded by a constant (expected: no; if yes: revolutionary).
- Where: findings/5 ("candidate weighting for future work"); findings/9
  (open leads ranked, #2).
- Class: COMPUTABLE (statistical test; a bound would then need proof).
- Unlocks: closes or spectacularly reopens the height-function route.

### 13. Structured-stream realization sweep at greater depth
- Statement: run the Realization Computer on more automatic/substitutive
  climbing streams and deeper digit windows (currently: TM and friends show no
  integer-like realization through 2^220).
- Where: findings/32 (R25, "C2 operational for arbitrary explicit streams");
  REPORT R201-205.
- Class: COMPUTABLE.
- Unlocks: the evidence base and counterexample-search for item 24.

### 14. Erdos-cousin scan: ternary digit 2 in 2^n beyond n = 4000
- Statement: extend the verification that 2^n contains ternary digit 2 for all
  n > 8 (2-free exponents exactly {0,2,8} so far).
- Where: findings/36 (R214-215, "the family's purest open problem").
- Class: COMPUTABLE (range extension only; the full problem is item 37).
- Unlocks: numerical territory on the closest sibling problem.

### 15. B mod M obstruction closed forms (Track A3)
- Statement: in Z/(2^S - 3^K), compute closed forms for B(shape) mod M, test
  equidistribution over falling shapes, hunt small-prime obstructions (the one
  observed skew was a size artifact -- question remains open, low priority).
- Where: PLAN Track A3 (P ~ 10%); findings/32 (R17).
- Class: COMPUTABLE.
- Unlocks: any structure found feeds the cycle-exclusion machinery (item 7).

### 16. Blocker-lattice conditional per-prime congruences
- Statement: develop the parked conditional per-prime congruence lattice into
  rigorous (conditional) exclusion instruments; note the unconditional blocking-
  primes mechanism itself was REFUTED (all primes reachable) -- only the
  conditional variant is still live.
- Where: findings/35 ("parked as rigorous future work"); REPORT 4f (blocking
  primes refuted, closed).
- Class: COMPUTABLE (exact automaton reachability at chosen moduli).
- Unlocks: conditional congruence certificates for orbit segments.

### 17. Lean formalization of the elementary theorems
- Statement: submit Theorems 1-3, 5, 9, 9', 10, 17, 20 (all formalizable) to
  the ccchallenge.org Lean pipeline.
- Where: PLAN C1; COMMUNITY.md (371 papers catalogued, 1 formalized).
- Class: COMPUTABLE (formalization labor, no new mathematics).
- Unlocks: machine-checked base layer; community adoption of the framework.

---

## LOCAL-ANALYTIC

### 18. Prove the interior cascade damping ratio is uniformly < 1
- Statement: the mid-cascade CV damping per digit level is measured <= 0.86
  across k = 13/15/17/20; prove it uniformly < 1 -- local statement (i) of the
  alpha -> 1 program.
- Where: NOTE.md Prop 23 ("two remaining ANALYTIC statements... (i)");
  REPORT R577-585.
- Class: LOCAL-ANALYTIC (explicit operator; Lemmas A/B of Thm 15 proven;
  1-Lipschitz minima; the attack surface is fully specified).
- Unlocks: half of alpha -> 1; with item 19 gives q -> 1, hence gamma -> 1 by
  Thm 19, hence pi(x) >= x^(1-eps) for every eps.

### 19. Prove CV_1 is bounded (fine-end saturation)
- Statement: the finest-level triple-CV satisfies CV_1(k) = 0.5136 -
  0.337*(0.910)^k (residual 1e-6, nine depths) -- prove the limit exists /
  the source is bounded -- local statement (ii) of the alpha -> 1 program.
- Where: NOTE.md Prop 23; findings/41; REPORT R586-594 ("condition (ii)
  closed numerically").
- Class: LOCAL-ANALYTIC.
- Unlocks: the other half of item 18's implication chain.

### 20. Prove the per-level injection amplitude decays
- Statement: linear damping is proven overwhelming (spectrum exactly {1} union
  {flat 1/4}, Thm 21); the g-field survives only by per-level injection through
  the min-term -- prove the injection decays and alpha -> 1 follows.
- Where: NOTE.md Thm 21 ("the alpha->1 question becomes: prove the injection
  amplitude decays -- damping is already proven overwhelming").
- Class: LOCAL-ANALYTIC (alternative, possibly cleaner route to items 18-19;
  Prop 18 linearizes the min-term as a single O(CV) correction).
- Unlocks: alpha -> 1 in one statement instead of two.

### 21. Derive the lattice coefficients (a,c) from the weight structure
- Statement: derive (a,c) of the CV recurrence CV_P = a*CV_(P-1) + c*CV_(P+1)
  analytically from the exact lattice identity (they are correlation-weighted
  masses), and show a - c stays bounded away from 0 past the k ~ 17 crossing.
- Where: NOTE.md Thm 16 ("Open to complete the proof: derive (a,c)... and show
  a-c stays bounded away from 0"); REPORT final synthesis ("two lemma-shaped
  steps").
- Class: LOCAL-ANALYTIC (exact identity verified to 2.8e-4; but note R206-210:
  the channels are uncorrelated, so no sign-coherence shortcut -- the derivation
  must go through the full spatial correlation structure).
- Unlocks: theta_inf < 1 from first principles; completes Theorem 16.

### 22. Prove the min-mean gap law (Prop 18) analytically
- Statement: E[1 - min/mean] = c1*CV + c2*CV^2 with c1 -> ~1.19 and c2 bounded
  in [-1, -0.5] -- currently MEASURED at 2 depths, 9 3-adic levels; prove it.
- Where: NOTE.md Prop 18; findings/41 (R409-420, "the linearization loop now
  has all its ingredients measured").
- Class: LOCAL-ANALYTIC (order statistics of a nearly-homogenized triple).
- Unlocks: rigorous linearization of the K-L min-nonlinearity -- the load-bearing
  lemma inside items 18, 20 and 21.

### 23. The Homogenization Conjecture / uniform spectral gap (q -> 1)
- Statement: prove a uniform spectral gap (subleading ratio bounded < 1) for the
  within-triple difference operator -- equivalently CV -> 0, q -> 1, alpha -> 1
  (Conjecture T) -- the umbrella assembled from items 18-22.
- Where: NOTE.md Thm 15c + Homogenization Conjecture + Conjecture T; findings/35
  (three instruments agree ~0.84); findings/40 ("OPEN CORE (new form): prove
  alpha_k -> 1"); REPORT R121-130, R445-472 ("the open core of the entire
  density program reduces to one statement").
- Class: LOCAL-ANALYTIC by the corpus's own assessment (every ingredient is
  lemma-shaped and measured), but the hardest item in this class; if the
  correlation structure resists, it escalates to HARD.
- Unlocks: pi(x) >= x^(1-eps) for every eps -- the full density version of the
  3x+1 problem, the K-L method's original hope.

### 24. The Mahler/Christol rung: no automatic climbing stream realizes a
### positive integer
- Statement: prove that the 2-adic realization of any automatic (or
  substitutive) climbing index stream is never a positive integer; the easy
  route (transporting automaticity through the transducer) is CLOSED -- the
  proof must work on the parity side via Mahler functional equations.
- Where: PLAN C2c ("concrete, novel, attackable"); REPORT R201-205 (rung-1
  verdict) and final synthesis ("first rung, buildable now").
- Class: LOCAL-ANALYTIC (Christol/Cobham/Mahler technology exists; novel
  application).
- Unlocks: first unconditional per-orbit rigidity class beyond eventually-
  periodic (Thm 9'); kills every automatic conspiracy; a publishable theorem.

### 25. Theorem 9 extension: periodic skeleton with bounded defects
- Statement: extend the shadowing bound from eventually-periodic streams to
  streams that are periodic up to finitely many (bounded) defects -- "finite
  unions, should work".
- Where: PLAN C2b.
- Class: LOCAL-ANALYTIC (same 2-adic separation mechanism, finite unions).
- Unlocks: widens the unconditionally-dead escape-mode class; stepping stone
  between Thm 9' and item 24.

### 26. Base-3/2 grammar completion (AFS carry conventions)
- Statement: the insertion law (odd T-steps insert a digit '1' in greedy
  base-3/2) is exact on orbit values but only ~33% on arbitrary odds, and even
  steps are not single edits -- complete the grammar under the Akiyama-Frougny-
  Sakarovitch carry conventions.
- Where: findings/34 (Base-3/2 Insertion Law); findings/35 (R131-134 PARTIAL,
  "full grammar deferred"); REPORT (representations worth future work: rational
  base 3/2, balanced ternary, Ostrowski w.r.t. log2(3)).
- Class: LOCAL-ANALYTIC (existing AFS literature; bookkeeping-heavy).
- Unlocks: a third exact language for the dynamics; candidate new invariants.

### 27. Derive the 3n+1 fingerprint analytically
- Statement: the digit-energy cascades differ between maps -- 3n+1 =
  (0.814, 0.059, 0.030) vs 3n-1 = (0.696, 0.178, 0.026) -- the first measured
  quantities that separate the maps at measure level; derive them from the
  mirror r -> -r acting on 3-adic digits.
- Where: REPORT R516-521 ("the place where a proof must live").
- Class: LOCAL-ANALYTIC (finite 3-adic computation + closed-form roulette).
- Unlocks: a theorem-shaped map-separating quantity -- the litmus-test
  requirement (must distinguish 3n+1 from 3n-1) finally met by an analytic
  object.

### 28. Dangling: restate the sign-coherence lemma post-refutation
- Statement: the final synthesis still lists "branch-channel sign-coherence" as
  a needed lemma-shaped step, but R206-210 REFUTED sign-coherence (channels
  uncorrelated, incoherence factor 0.90, branch:transport 4.2:1); the needed
  statement must be reformulated in the mean-field + fluctuation-coupling
  blueprint of R211-213 -- currently nobody has written down the corrected
  target lemma.
- Where: REPORT R206-213 vs REPORT final synthesis (inconsistency); findings/36.
- Class: LOCAL-ANALYTIC (a formulation task with measured constraints:
  uncorrelated channels, 0.90, 4.2 ratio).
- Unlocks: the correct target statement for item 21; removes a live
  inconsistency in the corpus's own proof-program description.

---

## HARD

### 29. H(eps): per-orbit quarter-fair mixing (Theorem C's hypothesis)
- Statement: prove that every orbit eventually staying above 2^71 satisfies
  freq(k_i >= t) <= (1+eps)*2^(1-t) and mean min(w_i,8) >= (1-eps)(2-2^-7) for
  some eps < 1/4 -- measured eps ~ 0.0002-0.01, required 0.25.
- Where: NOTE.md Thm 7; findings/24 ("the margin between measured and PROVEN is
  the conjecture"); REPORT section 5.
- Class: HARD (per-orbit equidistribution of the x3 action on 2-adic residues;
  Furstenberg x2x3 genre).
- Unlocks: the ENTIRE Collatz conjecture -- both divergence and cycles -- via the
  already-proven Theorem C.

### 30. E-star / the per-orbit reload law (Theorem 4's condition)
- Statement: prove that along EVERY orbit limsup (1/m) * Sum(k_i*log2(3) - k_i
  - w_i) < 0 -- equivalently that no positive integer realizes an infinite
  net-climbing stream (the archetype realizes at n = -5).
- Where: NOTE.md Thm 4 + "honest status" ("it IS the open problem") + exact
  restatement 2; findings/15 (Conjecture E-star); findings/16.
- Class: HARD (same wall as item 29 in martingale form; the machine is proven
  maximally forgetful -- Thms 20/21/22 -- yet per-orbit control is untouched).
- Unlocks: the divergence half; with the martingale identity, optional stopping
  finishes no-divergence.

### 31. C3: window-unfairness bound with ANY explicit f(v)
- Statement: prove "no orbit at value v is >= 26% unfair for m consecutive
  macro-steps once m >= f(v)" for any explicit f -- needs effective
  equidistribution of 3^j mod 2^i along orbit sequences (Baker/Ridout/
  exponential-sum territory). Envelope measured: unfairness never seen past
  m ~ 40-50.
- Where: PLAN C3 ("the precise, minimal statement of the whole problem");
  findings/32 (R29 window-unfairness envelope).
- Class: HARD (but graded: ANY f gives quantitative divergence-rate bounds;
  f = O(log v) gives the full conjecture -- the only known partial-credit path
  through the wall).
- Unlocks: partial credit at every strength level, up to the full conjecture.

### 32. The cycle half at scale: divisibility beyond the Catalan gifts
- Statement: prove no net-falling shape beyond (1,1) satisfies
  (2^(K+W) - 3^K) | B(shape) -- the genuinely open region is K > 6.5e10 with
  m > 91 circuits, where linear-form bounds degrade with term count.
- Where: NOTE.md Thm 5 + exact restatement 1; findings/8 ("open frontier");
  findings/20 (Baker-throttled coincidences).
- Class: HARD (effective transcendence at scale).
- Unlocks: the cycle half of the conjecture.

### 33. Effective subspace theorem for {2,3}-unit sums (two-prime case)
- Statement: prove an effective version of the subspace theorem for S-unit sums
  over {2,3} -- the first rung of item 32 named by the corpus.
- Where: REPORT final synthesis ("what crossing the wall requires... for
  cycles").
- Class: HARD (frontier Diophantine approximation; the two-prime special case
  is the least-hard entry point).
- Unlocks: converts the Baker throttle into finite computations -- could reduce
  the cycle half to a census run.

### 34. Terras CST conjecture (coefficient stopping = stopping)
- Statement: prove tau(n) = t(n) -- that the coefficient stopping time equals
  the stopping time -- the community's flagged clean open subproblem.
- Where: COMMUNITY.md ("clean open subproblem, community target").
- Class: HARD (open since 1976; but strictly weaker than the conjecture).
- Unlocks: community-recognized partial result; would sharpen every survivor-
  class computation in the corpus (F7, Thm 10).

### 35. Dynamic fuel accounting / a non-representational Lyapunov function
- Statement: find a quantity updated along the orbit (not read off the
  representation) that upper-bounds future ternlen streaks -- static positional
  scores are PROVEN insufficient (bounded-window impossibility, F11), scalar
  and phase-parity symbolic weights are proven impossible.
- Where: findings/6 (positional-fuel lead); findings/11 ("what remains open:
  dynamic fuel accounting"); findings/9 (Lyapunov equivalence); findings/34
  (scalar/phase-parity closures); REPORT R195 (the perfect Lyapunov function =
  steps-to-1, existence-as-finite = the conjecture).
- Class: HARD (the Lyapunov equivalence makes any success tantamount to the
  conjecture; every representational shortcut is now provably closed).
- Unlocks: the conjecture, via a certificate object.

### 36. The conjecture-equivalent reformulations (no independent handle yet)
- Statement: three exact restatements each awaiting any tool: (a) an
  unconditional peak cap f(n) (equivalent to the divergence half); (b) pointwise
  boundedness of the wavefront laggard slope c*(n) < infinity; (c) "measure zero
  is not empty" -- annihilating the survivor residue classes that thin as
  t^(-3/2)*2^(-t/20) but never empty.
- Where: (a) REPORT 4f peak-height answer, part 5; (b) findings/39 R354-358
  ("proving any uniform lower bound on the laggard-edge slope IS the
  conjecture"); (c) findings/7, findings/9 (the one irreducible core).
- Class: HARD (each IS the conjecture; listed so no thread is lost).
- Unlocks: everything; these are the wall itself in three coordinates.

### 37. The sibling family: Erdos 2^n-ternary, Mahler Z-numbers, Furstenberg x2x3
- Statement: the corpus identifies Collatz as one of four famous open problems
  sharing one core -- the arithmetic incommensurability of bases 2 and 3
  (guarded by Cobham's theorem: no finite-state dictionary exists).
- Where: REPORT R193 (conversion thesis), R214-215, final synthesis.
- Class: HARD (each sibling is itself a famous open problem).
- Unlocks: any breakthrough on any sibling transfers to the whole family;
  progress on item 24 or 33 is progress here.

---

## Summary count
- COMPUTABLE: 17 items (1-17)
- LOCAL-ANALYTIC: 11 items (18-28)
- HARD: 9 items (29-37)
- Total: 37 unsolved threads.

The corpus's own proof program, restated: items 18+19 (or 20) + 21/22 close
item 23 (q -> 1), Theorem 19 transfers it to gamma -> 1 (done), giving
pi(x) >= x^(1-eps); item 1 referees the whole chain. The full conjecture then
still requires item 29 (or 30/31) for divergence and item 32 (or 33) for
cycles.
