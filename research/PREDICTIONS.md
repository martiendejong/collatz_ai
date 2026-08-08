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

13. (frozen 2026-08-06, Obs 502) gamma-ladder law (1-gamma) ~ C*(2/3)^(k/6):
    gamma(25) = 0.9378 +- 0.002; gamma(27) = 0.9457 +- 0.002;
    gamma crosses 0.95 at k = 28.3 +- 1. (DENSITY fork; CEILING predicts flattening.)

14. (frozen 2026-08-06, Obs 503) true-lambda* extrapolation (Script 307):
    gamma_true(25) = 0.9429 +- 0.004; lambda*(21) = 1.8955 +- 0.003
    (=> polished k=21 certificate at lambda~1.894 should be feasible, gamma ~ 0.922).
    Rate of extrapolated (1-gamma): 0.926/step (~2^(-1/9)); certified-points rate
    0.9347 (~(2/3)^(1/6)); discrepancy attributed to lower-bound bias of certificates.

15. (frozen 2026-08-08, Obs 522/U2) s(k) := minus slope of linear fit of rho(lam,k)
    over lam in {1.85, 1.87, 1.89, 1.91, 1.93}. Prediction: s(k) decreases
    monotonically toward s(inf) = -rho_lin'(2) = 0.1556 (algebraic, from the mass
    identity). Measured so far: s(9)=0.3320 ... s(16)=0.2929, ratio ~0.982/step.

16. (frozen 2026-08-08, Obs 524) TR(2.00, k=16->17) = 0.8215 +- 0.0045, where
    TR = inc_last(k+1)/inc_last(k) per Script 333 (log2-field, prefix increments,
    last layer). Stationarity claim: no upward creep in the clean instrument.

    SCORECARD #16: MISS (measured 0.83005). Reinterpretation in Obs 525:
    series was still converging; limit matches rate(gammabar)^2 at 4 decimals.

17. (frozen 2026-08-08, Obs 525) TR(2.00, 17->18) in [0.828, 0.838] (center 0.833):
    the series has converged near c(2) = 0.830 = rate(gammabar)^2; no further rise
    beyond 0.838. Same definitions as #16 (Script 333/334 recipe, full convergence).
