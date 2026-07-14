# F40 — Campaign VIII (Rounds 359–408): genealogy, excluded set, tempering law

## R359–366 Champion genealogy
- 30% of A006877 records = backward propagation ((4n−1)/3 "+2 chains", 2n)
- seeds bimodal: siblings (share 70–97% of orbit with prev champion) vs true
  discoveries (share 2–9%); ALL true discoveries merge at value 20 = turnstile corridor

## R367–373 Excluded-set anatomy
- P(excluded|k) monotone in trailing-ones k at every t
- mod-64 concentration: 2.1x (median) → 75x (p99); top class always 63=111111₂, then 31,27,47

## R374–380 The γ–λ law
- (1−γ_k)/CV_res(k) = 0.6991/0.6973/0.6971/0.6985 — constant to 3 decimals
- density and spectral programs unified; prediction k=21: γ ≈ 0.918

## R381–408 The tempering law
- g ≈ −log(roulette) (r=−0.947); full linear model R²=0.94; rung1-spring r=+0.52
- eigvec = roulette^α: R²=0.9927/0.9949/0.9963/0.9973; α=0.8024/0.8291/0.8509/0.8682
- 1−α_k = CV_res(k) numerically; grand chain: γ_k = 1 − 0.698·(1−α_k)
- OPEN CORE (new form): prove α_k → 1 (approach ratio ~0.87/2digits, slight creep = θ∞ caveat)

Scripts: 48_champion_genealogy.py, 49_excluded_anatomy.py, 50_g_regression.py (+ inline power-law fit)
