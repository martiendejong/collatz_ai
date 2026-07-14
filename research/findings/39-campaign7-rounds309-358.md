# F39 — Campaign VII (Rounds 309–358): factorization law & backward wavefront

## R309–315 Exact roulette theory
- closed form mod 9: π(1,2,4,5,7,8) = (8,16,11,4,2,22)/63 ✓ (numeric 4-digit match)
- bridge triangle mod 9/27/81: theory↔orbits 0.997/0.997/0.992; theory↔eigvec 0.996/0.994/0.992

## R316–322 The cascade is the roulette
- roulette digit-energy: 88.4% digit 0, ratios 0.082/0.514/0.149 = eigvec spectrum incl. period-2
- residual CV 0.197 vs raw 0.90 → spectral content isolated

## R323–333 Factorization law + k=19 referee
- eigvec ≈ roulette × (1 + λ_k·g); g universal (cross-k corr 0.9999)
- λ: CV 0.197→0.171→0.150 (ratio ~0.874/2digits); k=19 PREDICTED 0.131, MEASURED 0.1333 ✓, shape corr 0.9991
- caveat: ratio creeps 0.871/0.876/0.887 — θ∞ question relocated

## R334–340 g identification
- NOT chain modes (|λ2|=0.003, R²=0.20); corr(g, first-preimage-is-spring)=+0.52

## R341–353 Backward wavefront (user's program)
- exact d(n) all n<2^24; laggard records = OEIS A006877 exactly ✓
- laggards are ladder-riders: mean k=3.27 vs 2.00; champions k=6..11 (26623=family 13, k=11)
- wedge edges: slope 1 (doubling) vs c*=d/log2n rising 12.6→18.4 (LW stochastic const ~29)
- figure: figures/wavefront_shape.png; same shape in graph/collatz_backward.graphml (depth attr)

## R354–358 Synthesis
- wavefront reformulation exact: bounded c* pointwise ⟺ conjecture
- laggard pattern locates extremal candidates: heavy ladders in small families

Scripts: 41_roulette_theory.py … 47_wavefront_shape.py
