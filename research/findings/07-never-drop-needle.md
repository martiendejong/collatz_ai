# F7 — The Never-Drop Needle: exponential thinning, survivors are trailing-ones rich

**Status: PROVEN (k=1 drop; DP counts exact) + MEASURED (survivor census).** Figure: `figures/e7_needle.png`

## The one-step theorem
Index k=1 (n ≡ 1 mod 4) forces a drop below self: n = 4b+1 → 3n+1 = 12b+4 → /4 → 3b+1 < n. Verified: all n ≡ 1 mod 4 up to 20,000 drop. **Minimal element of any cycle/divergent spine must have k ≥ 2 (n ≡ 3 mod 4)** — sole exception: 1.

## The lifelong constraint
A number never drops below itself iff at every prefix of its T-orbit, odd-step count K_j ≥ j / log₂3 = **0.63093·j**. Not a single index threshold — an infinite prefix inequality (one low stretch, ever, disqualifies).

## Exact survivor counts (Terras bijection: residues mod 2^t ↔ parity vectors)
| t | surviving fraction |
|---|--------------------|
| 20 | 2.61e−2 |
| 40 | 5.82e−3 |
| 80 | 6.64e−4 |
| 120 | 9.88e−5 |

Thinning rate: **2^(−0.050·t)** (KL divergence between coin 1/2 and threshold 0.631). Vanishing forever, empty never — at t=120 there remain ~10^32 surviving residue classes, each with infinitely many integers. **Measure zero ≠ empty is the entire remaining problem.**

## New structural finding: survivors are trailing-ones rich
Census of all odd n < 2^22 surviving 40 T-steps (23,924 survivors):
> mean trailing-ones of survivors = **4.94** vs **2.00** for all odd numbers.

Never-drop candidates are 2.5× enriched in the family/pair framework's own coordinate — the needle's survivors are exactly the high-sequence-index population, quantitatively confirming that the (a,k) representation is the natural habitat of the hard cases.

Related: [[05-heights-and-drift]], [[06-ternary-streaks]], [[09-verdicts-and-open-core]]
