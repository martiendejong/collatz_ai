# Concept-mail aan Jeffrey C. Lagarias (en cc Ilia Krasikov)

**Status: CONCEPT — versturen alleen door Martien, na eigen review.**
**Adres:** lagarias@umich.edu (University of Michigan; verifieer actueel adres op
zijn homepage vóór verzending). Krasikov: Brunel University London.

---

**Subject: Your L_k^NT linear programs solved for k = 12..18: pi_1(x) > x^0.902**

Dear Professor Lagarias,

In your paper with Ilia Krasikov, "Bounds for the 3x+1 Problem using
Difference Inequalities" (Acta Arith. 109 (2003), 237-258), the linear
program family L_k^NT(lambda) was solved up to k = 11, yielding
pi_1(x) > x^0.84, and the closing section expressed the hope that the
family might eventually be pushed further.

I write to report that I have solved the family for k = 12 through 18 on
modern hardware. The resulting exponents are:

    k = 12: 0.8520   k = 15: 0.8801
    k = 13: 0.8620   k = 16: 0.8875
    k = 14: 0.8715   k = 17: 0.8950
                     k = 18: 0.9020

each backed by an explicit feasible solution of L_k^NT(lambda) (certificate
vector with verified margin min F(v)/v >= 1.0004, coefficients computed at
60 digits and rounded downward), so that by your Theorem 2.2 each yields
pi_1(x) > x^{log2 lambda} unconditionally. The implementation was calibrated
against the published values at k = 2, 9, 11, and the index algebra verified
against the formulas of Proposition 2.1 directly.

Two questions, if I may:

1. Are you aware of any computation of this family beyond k = 11 in the
   intervening years? I have found none in the literature, but you would
   know best.

2. The sequence gamma(k) fits 1 - gamma(k) ~ C q^k with q about 0.93,
   consistent with gamma(k) -> 1 (the hope stated at the end of your
   paper). In the course of this work I also found what appears to be
   exact structure in the certificate vectors: the backbone m -> 4m is a
   single N-cycle, the feed maps contract 3-adic agreement by exactly one
   digit, and the total feed share of the Perron flow equals 1 - lambda^-2
   identically (the backbone being a permutation). Might this be of
   interest, or is it known?

All code, certificate vectors, and a full research log are public at:
https://github.com/martiendejong/collatz_ai
(directory research/, scripts 163-166; certificates as .npy files;
draft write-up in research/draft-arxiv-note.md).

I am an independent researcher (no institutional affiliation); the
computations were carried out with AI assistance and every step is
reproducible from the repository. If the result is new, I would be
grateful for your advice on an appropriate venue — and, should you be
willing, an arXiv endorsement for math.NT.

With admiration for your work on this problem,

Martien de Jong
[woonplaats], The Netherlands
[e-mailadres]

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
