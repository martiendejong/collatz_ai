# F42 — Campaign XIV (Rounds 701–835): CST, the master constant, the drift arc

## R701–710 (a−c) flow
- negative and moving away from 0 under every fit convention; flow fit -> −0.68 (quick-fit conv.)
- caveat: absolute (a,c) differ per windowing convention; direction robust

## R711–745 CST arc (Props 28, Thm 29)
- CST verified tau<=24 (81,119 classes; only n=1); margins governed by (5,8) = Pythagorean LIMMA 256/243
- direct basis: tau=sigma all odd n<=10^6 (0 violations, max tau 176)
- reduction: violation needs n <= b/(2^j−3^u); B-bound b/2^j < u/3 PROVED (paper) -> CST tau<=4700
- corrected: 2^65/3^41 = 1.011529; first convergent escape (15601,24727), threshold 2.86e8
- naming: 256/243 = limma; comma proper = 3^12/2^19

## R746–755 B-growth (Lemma 30)
- B(t) ~ 1 + t/8; extremal classes = REPEATED limma words (5,8),(10,16),(15,24)

## R756–765 Master constant (Prop 31)
- delta = log2(16/9) = 2log2(4/3) = 0.830075: flow rate = delta (0.01%), saturation = sqrt(delta) (0.12%)
- edge rate = 1/ln(4/3) closed form (paper upgrade); theta_inf/CV1_inf: no closed form yet

## R766–785 Drift skeleton (Prop 32)
- BALANCE IDENTITY PROVED: up-flow = W0*3/4 = 3/16 = wbar*1/4 = down-flow (third exact balance from 2^alpha=3)
- kappa(P) = 0.908→0.841 (P=2..7, k13) converging into theta range (~0.849)
- kappa uniformity k=15/17/19: background job launched (b9cdngd99; k19 memmap slow, still pending)

## R786–805 Gaussian reduction FAILS
- MC kappa ~0.98-1.00 vs measured 0.84-0.91, even with mean-reversion rho_LI = −0.085..−0.184
- kappa is higher-order; lesson: exact identity route only

## R806–815 Clipping decomposition (Prop 33)
- kappa^2 = Var(Dbar)/V + 2cov/V + Var(R)/V exact; ALL attenuation in cov (−0.17→−0.26)
- R antisymmetric ±0.026; effective Dmin = (1−lambda)Dbar; lambda = 0.085→0.137

## R816–825 Switch-resolved (Prop 34)
- asymmetry NOT switch-exclusive (P(sw)~0.66 flat) => selection-weighted mean reversion
- LP interpretation: argmin member = binding constraint (complementary slackness rigidity)
- proof route: LP duality / sensitivity

## R826–835 Directional low-pass (Prop 34b, proved core)
- min(x+c·1) = min(x)+c: coarse modes slope EXACTLY 1; intra-triple modes contracted (1−lambda)
- drift c>a = one-sided low-pass: down-channel starves fine modes, up-channel free
- remaining: contraction uniform < 1 (via Prop 23 nonzero saturation) + recurrence bookkeeping

## Also this campaign
- FOUR PAPERS completed (density 0.9146, tempering_law, cst_comma, forgetful_machine) + 3 math upgrades from writing
- PREDICTIONS.md (13 pre-registered), RUNBOOK_k21.md
- gamma->1 chain final state: Lemma24(P) + balance(P) + [min low-pass contraction uniform<1] + Thm19(P)
