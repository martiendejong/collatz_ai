# F25 — Records Scan to 2×10⁸: Prediction Falsified, Tail Rate is Scale-Dependent

**Status: MEASURED (vectorized scan, n < 2×10⁸; deeper scan to 2×10⁹ running).**

## True record chain (post-filtered, σ strictly increasing in n)
27:59 · 703:81 · 10087:105 · 35655:135 · 270271:164 · 362343:165 · 381727:173 · 626331:176 ·
1027431:183 · 1126015:224 · 8088063:246 · 13421671:287 · 20638335:292 · 26716671:298 ·
**56924955:308** · **63728127:376**

(8088063 and the famous 63728127 — the celebrated deep-delay number — recovered blind. ✓)

## The falsification
F23's linear law (slope 9.5, from records ≤ 1.1×10⁶) predicted first σ ≥ 300 near n ≈ 4×10⁹.
**Actual: n = 56,924,955 — seventy times earlier.** And 63728127 runs at σ/log₂n = 14.5,
above even the naive tail cap 1/0.080 = 12.5.

## Resolution: the tail exponent bends toward the KL floor
The stopping-time tail rate is scale-dependent:
- mid-range (t ∈ [20,140], measured): 0.080 bits/step
- implied by the deep records (σ = 376 among ~1.6×10⁷ candidates): ≈ 0.064 bits/step
- theoretical endpoint floor (F7 KL rate): 0.050 bits/step

The rate decays toward the floor as t grows — deep survival is asymptotically *cheaper per step*
than mid-range survival (large-deviation prefactors; paths hugging the 0.631 line dominate).
Consequently record growth accelerates relative to any mid-range linear fit; the safe asymptotic
envelope is σ_max(N) ≤ log₂N / 0.050 = 20·log₂N.

## Lesson (methodological, kept for honesty)
A 10-point fit in the pre-asymptotic regime produced a falsifiable prediction that failed within
one session. The failure was more informative than the fit: it exposed the scale-dependence of the
survival price — a real feature of the large-deviation structure that mid-range statistics hide.

Related: [[22-period4-and-records]], [[07-never-drop-needle]], [[23-synthesis-consistency-web]]
