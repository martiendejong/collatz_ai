# F38 — Campaign VI (Rounds 259–308)

## R259–265 Cascade spectrum (cert k13/15/17)
- ~89% of Var[log C] in finest digit (mod-3/9 branch structure)
- energy decay ratio ≈ 0.20/digit, stable across k (0.199/0.200/0.222), period-2 modulation
- total Var grows slowly: 0.663/0.713/0.757

## R266–271 Spine anatomy
- spine ≠ level minimum (depth-9: spine 3077 vs min 33)
- spine = consecutive forward trajectory: the 27-highway (911,1367,2051,3077→9232 peak)
- corridor entropy ≈ 0.3 bits/layer → effective branching 1.23

## R272–277 Multifractal spectrum
- D1≈0.92, D2≈0.81, D4≈0.67, D∞≈0.52; D∞ = spine exponent 1/2 (matches R258)

## R278–283 Caste Markov
- mod 3: memoryless (1/3, 2/3) exact; caste = w-parity
- mod 9: stationary 8:.356 2:.242 4:.172 1:.124 5:.072 7:.034; ord9(2)=6 roulette
  next ≡ (3m+1)·5^w; |λ2|=0.085 (memory dies in 0.28 steps)

## R284–289 Survivors + first bridge
- S(k) ~ 1.834^k (density 0.917^k), period-3 ratio cycle 2.00/1.81/1.70; S(16)=2114
- eigenvector mod-9 means track forward visit shares within 10%

## R290–295 Graph census
- tower-law in-degree formula: 200/200 exact
- ternary alternators (9^j−1)/8 odd j: in-degree ranks 3/6/7 of 79,545
- caste duality: forward-heavy 8 is backward-poor (0.68) vs 4 (2.65)

## R296–301 Controls
- 3n−1: basins ≈ thirds (.328/.324/.349); same turnstile funnels per cycle (43+11=98%; 19; 163)
- 5n+1: gates (2^k−1)/5, 4|k; every 5th gate multiple of 5 = shut
- general: cn+1 gate castes cycle mod c, period c

## R302–308 Eigen-visit bridge + synthesis
- mod-27 comparison across 9 classes: r=0.985 (log-log 0.989)
- LP certificate fine structure = forward visit measure (time-reversal duality)
- k=20 referee still running (~iter 110/150)

Scripts: 33_cascade_spectrum.py … 40_bridge_mod27.py
