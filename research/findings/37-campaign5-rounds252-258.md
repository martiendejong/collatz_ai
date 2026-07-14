# F37 — Campaign V (Rounds 252–258): doorway anatomy & visit measure

## R252–253: Gateway census
- Spring gates (y≡0 mod 3: 21, 1365, …) receive EXACTLY zero traffic (200k orbits) ✓
- Gate 5 takes 93.8%; caste ≡2 mod 3 gates beat ≡1 gates (w=1 shrink-move → denser basin)
- Door-w is caste-quantized (odd w for g≡2, even for g≡1) with spring-blocked rungs:
  gate 5: w=1→3 SHUT, w=3→13 (0.478), w=5→53 (0.456), w=7→213 SHUT, w=9→853

## R254: Turnstiles
- Penultimate odd ∈ {13, 53}: 93.1%. Antepenultimate ∈ {17, 35}: 91.0%.
- Terminal corridor: {17,35} → {13,53} → 5 → 1.

## R255: Funnel law
- Top-1 layer mass plateaus 0.37–0.40 through depth 14; top-2 decays ×0.95/layer.
- Persistent spine; corridor half-life ≈ 13 layers.

## R256: Eigenprofile spatial structure (cert_k13/k15)
- CV field locally white (lag-1 −0.02, ξ=1), ultrametrically correlated at coarse
  3-adic scales (r→+0.25 at top-digit lag, k=13; +0.16 at k=15 same lag).
- corr(CV, mean) = +0.35. Spectral problem sharpened: coarse-digit coupling is the content.

## R257: Record stratigraphy (12235060455, σ=547)
- Climb mean k = 2.99 (random climbs 2.46, own descent 2.01 = stationary ✓), max k=11.
- Ladder burns k=11→10→9 (one rung/odd step, w=1); ladders ≈ 9% of climb steps.

## R258: Visit measure
- Mean P(visit v) ∝ v^−1.01 exactly (complete windows, zeros incl.).
- Naive c≥50 fit gives −0.50 = truncation bias (documented trap).
- Multifractal: max/mean grows ~v^0.56; spine measure ~v^−0.44, bulk 1/v.

Scripts: 27_gateway_census.py, 28_turnstiles.py, 29_funnel_law.py,
30_cv_spatial.py, 31_record_stratigraphy.py, 32_visit_measure.py, 32b_visit_measure_clean.py
