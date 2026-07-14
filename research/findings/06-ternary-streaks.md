# F6 — Ternary Length: best Lyapunov raw material; streaks die geometrically, no cap found

**Status: MEASURED (20,000 orbits × 60 steps; champions exhaustive to 2^22).** Figure: `figures/e3_streaks.png`

## Why ternary length is the premier candidate quality
- **Provably bounded increase:** 3n+1 appends exactly one ternary digit; halving never lengthens. So ternlen rises by at most +1 per Syracuse step.
- Rises on only ~19% of steps; average drift −0.26 digits/step.
- Lives natively in the rewriting language (F2): each odd step writes one trit, each halving erases ~0.63 of one. Conjecture ⇔ "the eraser outpaces the pen on every orbit."

## Streak analysis (runs of consecutive non-decreases)
Continuation ratios P(streak ≥ j+1 | ≥ j) for j = 1..8:
**0.666, 0.682, 0.679, 0.677, 0.682, 0.687, 0.695, 0.694**

- Flat ≈ 0.68 → geometric death, **no structural cap detected**. (Slight upward creep suggests mild self-sustaining structure in survivors — consistent with F7's trailing-ones enrichment.)
- Champions to 2^22 (streaks 42–44): binary shows internal **alternation blocks** (…10101…), explained below; their *leading* trailing-ones are small (1–2), so the fuel is positioned deep, spent later.

## The alternation-fuel identity (exact)
> 3 · (01)^m ₂ = (1)^2m ₂  — i.e. 3·(4^m−1)/3 = 4^m−1: multiplication by 3 converts binary alternation into solid trailing-ones runs.

Alternation is **latent sequence-index fuel**: an orbit carrying 01-blocks in its low bits keeps reloading large k. Champions carry visible alternation blocks (e.g. 917161 = 11011111111010101001₂).

## Honest negative result
Global alternation density does NOT separate champions from average (0.4975 vs 0.4992, n=40,000). The fuel effect is positional and dynamic (low-order bits at the moment of use), not a static global score. A positional/weighted version remains open — this is the sharpest open experimental lead of the program.

Related: [[02-rewriting-theorem]], [[07-never-drop-needle]]
