# F36 — Campaign IV (Rounds 198–247)

## R198–200: the k=20 referee (launched)
Chunked float32 eigenvector at k=20 (1,162,261,467 classes, λ=1.885, 150 iterations, ~10 GB,
hours of runtime — the largest K–L computation ever attempted on this machine). Purpose: test
the lattice-model predictions ((a,c) ≈ (0.472, 0.521), θ ≈ 0.849) one depth beyond all fits.
Script: `scripts/25_k20_profile.py`; output: `k20_result.txt`.

## R201–205: Rung 1 verdict — the transducer destroys automaticity
Realization digits of automatic parity streams (Thue-Morse, period-doubling, paperfolding)
have subword complexity ≈ RANDOM (204–206 vs random's 203 at window 12; automatic level would
be ~20–48). Only periodic streams keep low-complexity (rational) realizations. Consequences:
(i) the easy Christol route (transport automaticity through the transducer) is CLOSED —
rung 1 must work on the parity side via Mahler functional equations; (ii) encouraging for the
deep goal: the transducer is complexity-EXPANSIVE, and every tested structured stream's
realization shows no trace of a positive integer's eventually-zero digits through depth 2^220.

## R206–210: sign-coherence hypothesis REFUTED
Transport and branch channels of the exact lattice identity are UNCORRELATED (corr −0.02/−0.05,
same-sign 48.7%/50.0% at k=13/17), incoherent-addition factor 0.90, branch:transport magnitude
4.2:1. The low effective mass is NOT sign-cancellation between channels — the lattice
coefficients emerge from the full spatial correlation structure, i.e. deriving them analytically
IS the spectral gap problem (no shortcut). The measured facts (uncorrelated channels, 0.90
incoherence, 4.2 ratio) constrain any future derivation.

## R211–213: the transmission mechanism identified
P(argmin unchanged under a mid-position offset) ≈ 0.33–0.38 — barely above random (1/3): the
min-selection RE-ROLLS almost completely under perturbations. Yet transmission is 0.8–0.9:
resolution — at deep k the triples are nearly homogenized, so the min tracks the triple MEAN
(c̄ ≈ mean·(1 − 1.36·CV)); differences transmit through the mean channel, robust to selection
switching. Future derivation blueprint: mean-field + fluctuation coupling, not selection
bookkeeping.

## R214–215: the Erdős cousin verified in range
2^n contains ternary digit 2 for all 8 < n ≤ 4000 (2-free exponents: exactly {0, 2, 8}) —
the family's purest open problem, consistent as far as tested.
