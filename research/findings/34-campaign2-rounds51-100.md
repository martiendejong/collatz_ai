# F34 — Campaign II (Rounds 51–100)

## New theorems
- **Theorem 11 (Residue-Blind Impossibility):** at every 2-adic depth j, class(−5) = 2^j−5 makes a
  consistent ×9/8 macro-step back into class(−5): no strict Lyapunov certificate on residue states
  exists at any finite depth. The negative cycles are proof obstructions for all congruence-state
  methods. (Verified j = 10–60.)
- **Theorem 12 (Min-Loss Identity):** 1 = λ^(−2) + (q/3)(λ^(α−2)+λ^(α−1)) exactly at every K–L
  Perron edge (verified 8 decimals, k = 5–11). q = 1 ⟺ λ = 2 ⟺ γ = 1: the method's gap to full
  density IS the min-loss.
- **Proposition 13 (Portrait of m★):** eight simultaneous unconditional properties of the minimal
  counterexample (in NOTE.md).

## The ceiling tension (open, decidable)
Two extrapolations of q_k disagree:
- Shanks on q-series (k = 5…19): q∞ ≈ 0.9927 → ceiling γ∞ ≈ 0.9757.
- CV-route: 1−q ≈ 1.36·CV(k) with intra-triple CV decaying geometrically (×0.82/two-depths,
  0.0368 → 0.0210 at k = 13→19, no floor visible) → **q → 1: no ceiling below γ = 1** (K–L's
  original x^(1−ε) hope). Projections under CV-model: γ(21) ≈ 0.927, γ(41) ≈ 0.988.
- Decisive test: k = 21 (3.5×10⁹ classes, ~16 GB — cloud task).

## New laws and measurements
- **Balanced-Ternary Anchor Law** (0/20,000): bal(E) = bal(a) | 0^(k−1) | T̄ — the anchor's native base.
- **Base-3/2 Insertion Law** (25/25): odd T-steps insert a single digit '1' in greedy 3/2-digits
  (×3/2 = shift there); even steps are not deletions. Third face of the two-language tension.
- **Homogenization mechanism:** q rises because Perron refinement-triples equalize with depth
  (near-equal triples 77% → 95%, k = 13→19).
- **Spectral gap trend:** edge-linearized ρ₂/ρ₁ = 0.64–0.87 at k = 4–9, creeping upward.

## Closures (exhaustive/proven negatives)
- Scalar symbolic weights: proven impossible (3-line proof). Phase-parity weights (6-dim, 262k
  candidates): best margin −1 — closed.
- Census: periods 8, 9, 10 (k,w ≤ 2) empty; 7-circuits (74.4M shapes) empty.
  Cycle list {1, −5, −17} complete through period 10 (small symbols) / circuits m ≤ 7.

## Infrastructure
- verify_certificates.py rewritten chunked/memmap (constant memory): k = 13/15/17 re-verified;
  k = 19 chunked run in progress. Records scan 8→20×10⁹ running.
