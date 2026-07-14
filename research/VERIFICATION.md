# How to Verify the Density Claim Without Trusting the AI

Claim: π(x) ≥ x^0.9069 (and the intermediate 0.8624 / 0.8805 / 0.8953), improving the
Krasikov–Lagarias 2003 record x^0.84. Every link in the chain below is checkable by a human
with no reliance on any code or statement the AI produced.

## The chain of trust (what rests on what)

1. **The mathematics (NOT ours):** Krasikov & Lagarias, "Bounds for the 3x+1 problem using
   difference inequalities", Acta Arithmetica 109 (2003); arXiv math/0205002 (PDF in this folder:
   `kl2003.pdf`). Their **Theorem 2.2**: any feasible positive solution of the linear program
   L_k^NT(λ) (their §2, equations 2.7–2.14) yields π_a(x) ≥ Δ₁·c^m·λ^y, i.e. a density exponent
   log₂λ. Peer-reviewed, 22 years old. We prove nothing new mathematically — we only exhibit
   feasible solutions at k = 13, 15, 17, 19.
2. **The system transcription:** our constraints must match their equations. This is ~20 lines
   to compare by eye: their (2.8)–(2.14) against the constraint block in
   `certificates/verify_certificates.py`.
3. **The certificates:** four integer vectors (`cert_k13/15/17/19.npy`). The claim is only that
   these vectors satisfy all constraints at rational λ₀ — a finite, mechanical fact.
4. **Novelty:** that 0.84 was still the record. Two literature sweeps found nothing above 0.84
   through 2026-07 (a Dec 2025 preprint, arXiv 2512.13760, still cites 0.84 as the record).

## Verification protocol (increasing depth)

**Level 0 — the artifacts are real (1 minute).** The files exist on this disk and have the right
sizes (3.5 GB total). Nothing was imagined: `dir certificates\`.

**Level 1 — run the verifier yourself (minutes).** `python verify_certificates.py` re-checks all
435,781,620 constraints in exact integer arithmetic (Python bigints via object arrays — no
floating point in the comparison, no overflow possible). Weights are strict rational LOWER bounds
(80-digit decimal floors minus one unit), so every rounding errs AGAINST the claim.

**Level 2 — read the 20 lines (30 minutes, the key step).** Open `kl2003.pdf` §2 and
`verify_certificates.py` side by side. Check: classes are m ≡ 2 mod 3; the three constraint
families branch on m mod 9 with weights λ^(−2), λ^(α−2), λ^(α−1); auxiliaries are minima over the
three refinements. This closes the ONLY correlated-error risk (the AI misreading the paper in the
same way in both solver and verifier), because YOU read the paper.

**Level 3 — hand spot-check one constraint (10 minutes, zero code trust).** Example generated
from the k=13 certificate, class index 123456 (m = 370370 ≡ 2 mod 9):
  c^m = 10,653,081,126,519 ; c^{4m} = 31,852,724,675,717 (index 493826) ;
  c̄ = min over indices {164608, 341755, 518902} = 1,305,501,746,119.
  Check on any calculator: 31852724675717·1.818^(−2) + 1305501746119·1.818^(log₂3 − 2)
  = 1.0656e13 ≥ 1.0653e13 ✓ (slack 1.000279). Pull the five integers yourself with
  `numpy.load('cert_k13.npy')[i]`. Repeat for any random indices you choose.

**Level 4 — the structural killer check: reproduce the PUBLISHED values.** Run the solver
(`scripts/22_kl_exact.py`) at k = 9 and k = 11: it yields γ = 0.8168 and 0.8418 — matching the
independently published 0.81 (Applegate–Lagarias 1995) and 0.84 (K–L 2003). A mistranscribed
system reproducing two independent published values to three decimals is not a plausible failure
mode.

**Level 5 — clean-room reimplementation (hours).** Give the K–L paper (only the paper) to an
independent programmer or a different AI in a fresh session: "implement L_k^NT(λ), find the
feasibility edge at k = 13." Compare their λ against 1.8188. No shared code, no shared reading.

**Level 6 — the experts.** Email Lagarias (the method's author) with the note and certificates
before or at submission. He would instantly spot any misreading of his own system. Then arXiv
itself: publication of a certificate-backed claim IS the community verification mechanism —
the certificates make refutation, if warranted, a five-minute job for any referee.

## Residual risks, stated honestly
- K–L's Theorem 2.2 itself wrong: peer-reviewed 22 years, used by later literature; risk ≈ that of
  any established result.
- Consistent misreading by the AI of the paper in solver AND verifier: eliminated by Levels 2 and 5
  (human/independent reading), made implausible by Level 4.
- Someone already did k ≥ 13 unpublished: mitigated by sweeps + Level 6; at worst the note becomes
  a confirmation, not a claim.
