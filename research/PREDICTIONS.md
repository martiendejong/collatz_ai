# Pre-registered predictions (frozen 2026-07-16, before any k=21 computation)

## k=21 K-L run (needs ~64GB; predictions frozen now)
1. gamma(21) ~ 0.918 (lambda* ~ 1.890); via (1-gamma) flow ratio ~0.885/2digits
2. alpha_21 = 0.887 +- 0.003 (tempering exponent, block mod 3^7)
3. CV_res(21) = 0.116 +- 0.004
4. theta(21) = 0.850 +- 0.001 (lattice fit, windowed convention of 25b)
5. (a,c)(21) ~ (0.465, 0.528) (Thm 16 corollary convention)
6. cascade digit-energy ratio ~ 0.20-0.22 per digit (fine half)
7. q(21) ~ 0.9755 ((1-q) ratio ~0.85/2digits from 0.97232 at k=20)

## Structural constants (any future depth)
8. (a-c) flow rate = delta = log2(16/9) = 0.830075 exactly
9. fine-end saturation rate = sqrt(delta) = 0.911084 exactly
10. CV_1(k) -> 0.5136 (saturation limit)
11. kappa (min-attenuation) at top-aligned deep scales -> theta_inf ~ 0.849-0.850
12. edge rate: (1-gamma)/( (1-q)/ln(4/3) ) -> 1 monotonically

## CST
13. no tau!=sigma violation will ever be found below the (15601, 24727)
    convergent zone (threshold 2.86e8) except n=1

## The gamma fork (added R2250, 2026-07-18)
- Both models predict gamma_21 = 0.919 +- 0.001 (cannot discriminate at k=21).
- CEILING model (gamma_inf = H(1/log2 3) = 0.950): gamma flattens, never crosses 0.950.
- DENSITY model (gamma_inf = 1): gamma crosses 0.950 around k ~ 27 +- 3.
- Discriminating experiment: exact-integer certifications at k = 25-30 (cloud scale).
