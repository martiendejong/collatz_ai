# F41 — Campaign IX (Rounds 409–508): the analytic skeleton

## R409–420 Min-mean gap law (Prop 18)
- gap = c1*CV + c2*CV^2; c1 in [1.19,1.45] -> ~1.19 as CV->0; c2 in [-1,-0.5]
- K-L nonlinearity = single O(CV) correction on roulette mean

## R421–432 Cycle census period 12 (Thm 17)
- complete integer cycle list through k=12 odd steps: {1,-1,-5,-17} (was: period 10)

## R433–444 Laggards to 2^25
- one new record: 31466382 = 2*15733191 (chain phase); c* = 17.75

## R445–472 Edge rate (Thm 19)
- dgamma/dq|_edge = 3.47614 EXACT from Min-Loss Identity
- (1-gamma)/(3.4761*(1-q)) = 0.824/0.847/0.873/0.908 -> 1 (k=13..19)
- 0.698 constant DERIVED as finite-k composite; open core = CV contraction only
- q values: 0.9519/0.9594/0.9655/0.9705; (1-q) ratio ~0.85/2digits

## R473–508 k=20 + synthesis
- k=20 first run lost (no checkpoint); rerun with checkpoints running
- proof program: (i) CV-cascade contraction [OPEN]; (ii) Thm 19 [done]; (iii) Thm C [conditional]

NOTE.md: +Thm 17, Prop 18, Thm 19, Conjecture T.
Scripts: 51_linearization.py, 52_census_p12.py, 53_laggards_2e25.py, 54_derive_0698.py
