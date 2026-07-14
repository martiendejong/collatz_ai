# F11 — Lookahead Quantified: predicting ascent costs bits linearly, determinism exponentially

**Status: MEASURED (30,000 orbits, 40–44 bit starts).**

## Experiment
Predict the 60-step ascent A(n) = log₂(max future value / n) from only the j lowest bits of n:

| j low bits | R² (variance explained) |
|---|---|
| 4 | 0.399 |
| 8 | 0.639 |
| 12 | 0.789 |
| 16 | 0.946 |
| 20 | 0.998 |

Mean ascent 2.71 bits. Full determination of a mere 60-step window needs ~20 bits of state.

## Reading
- Consistent with the lookahead regress (F1): the low bits ARE the future, but each deeper
  macro-step consumes them geometrically; r reloads of foresight cost ~2^r bits.
- The positional-fuel program (F6) is hereby bounded: any static weighting of low-order structure
  can control only an O(log j)-macro-step horizon. **A bounded-window positional fuel score cannot
  cap streaks globally** — closing open lead #1 negatively at the static level.
- What remains open: *dynamic* fuel accounting (a quantity updated along the orbit rather than
  read off the representation), which is no longer a "quality of the number" in the original sense.

Related: [[01-family-pair-coordinates]], [[06-ternary-streaks]]
