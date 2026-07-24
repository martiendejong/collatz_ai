# Concept-mail aan Jeffrey C. Lagarias (en cc Ilia Krasikov)

**Status: CONCEPT — versturen alleen door Martien, na eigen review.**
**Adres:** lagarias@umich.edu (University of Michigan; verifieer actueel adres op
zijn homepage vóór verzending). Krasikov: Brunel University London.

---

**Subject: Your L_k^NT linear programs solved through k = 20: pi_1(x) > x^0.9146**

Dear Professor Lagarias,

In your paper with Ilia Krasikov, "Bounds for the 3x+1 Problem using
Difference Inequalities" (Acta Arith. 109 (2003), 237-258), the linear
program family L_k^NT(lambda) was solved up to k = 11, yielding
pi_1(x) > x^0.84, and the closing section expressed the hope that the
family might eventually be pushed further.

I write to report that I have solved the family through k = 20 on a desktop
machine (3^19 = 1.16 billion constraints at the top level), obtaining

    k = 13: 0.8624   k = 17: 0.8953
    k = 15: 0.8805   k = 19: 0.9069
                     k = 20: 0.9146

each backed by a deposited explicit feasible solution of L_k^NT(lambda_0)
at rational lambda_0, re-verified in exact integer arithmetic with all
weights replaced by strict rational lower bounds (so the check errs
strictly against the claim), so that by your Theorem 2.2 each yields
pi_1(x) > x^{log2 lambda_0} unconditionally. Two independently written
implementations (different indexing and lambda parametrizations) agree on
the overlapping levels to 3-4 decimals, and the k = 12 certificate
(gamma = 213/250 = 0.852, already above 0.84) has additionally been
machine-verified in Lean 4 via exact integer power inequalities. The
implementations were calibrated against the published values at
k = 2, 9, 11, and the index algebra verified against Proposition 2.1
directly.

One question only: are you aware of any computation of this family beyond
k = 11 in the intervening years? I have found none in the literature
(a December 2025 preprint, arXiv:2512.13760, still cites x^0.84 as the
record), but you would know best.

Everything needed to check the claim without trusting me is public at
https://github.com/martiendejong/collatz_ai : the deposited certificate
vectors, a standalone numpy-only exact-integer verifier
(research/certificates/verify_certificates.py), the Lean 4 project
(research/lean/CollatzCert), and a three-page write-up
(research/NOTE_DENSITY.tex).

I am an independent researcher, assisted by Jengo, my AI research
assistant; every step is reproducible from the repository without trusting
either of us. Should you be willing, I would be grateful for an arXiv
endorsement for math.NT.

With admiration for your work on this problem,

Martien de Jong
Nijeveen, The Netherlands
martiendejong2008@gmail.com

---

## Toelichting bij het concept (niet meesturen)

- **Toon**: feitelijk, kort, geen grote claims — het woord "Collatz opgelost"
  komt er niet in voor. De twee vragen geven hem een makkelijke reden om te
  antwoorden (expertise-vraag, geen beoordelingsverzoek).
- **De structuurvraag (2)** is bewust compact: drie exacte feiten, geen
  theorie-verhaal. Als hij bijt, volgt de rest vanzelf.
- **AI-vermelding**: transparant en kort; verzwijgen is riskanter dan melden.
- **Endorsement-verzoek**: aan het eind, voorwaardelijk geformuleerd, geen druk.
- **Vóór verzending nodig van Martien**: woonplaats + e-mail invullen, actueel
  adres van Lagarias verifiëren, en eigen slotlezing.
