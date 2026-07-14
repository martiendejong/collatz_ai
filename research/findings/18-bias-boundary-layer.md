# F18 — The Borrow Bias is a Boundary Layer: 0.8/L, fully explained, no bulk structure

**Status: MEASURED + mechanism identified. Open lead #2 (round 3) closed negatively.**

## Size scaling
Deviation of P(0→1) from ⅓ on real cascades, by input size:

| bits | dev | dev × trit-length |
|---|---|---|
| 24 | +.0351 | 0.84 |
| 48 | +.0167 | 0.80 |
| 96 | +.0081 | 0.78 |
| 128 | +.0063 | 0.81 |

**dev ≈ 0.8 / L** — pure finite-size boundary effect, vanishing in the bulk limit.

## Localization (64-bit inputs, by trit distance from LSB)
- Position 0: **fully deterministic** (P=0) — the last-trit theorem: 3n+1 ≡ 1 (mod 3) always,
  and the final borrow equals x mod 2 = 0 for even x.
- Positions 1–5: strongly skewed (e.g. pos 1: P01=.46, P11=.91) — the LSB layer where the
  3n+1 structure and halving parity interact deterministically.
- Positions ≥ 6: **converge to (⅓, ⅔) exactly** (.333/.667 at pos 7–20).
- MSB end: mild skew (leading trit ∈ {1,2}, initial borrow = 0).

## Conclusion
The 2% deviation seen in F15 is the O(1)-digit boundary layer averaged over O(L) digits.
The bulk borrow process is perfect (⅓, ⅔) Markov to measurement precision. There is **no persistent
exploitable bias**: the equidistribution wall is reconfirmed at digit-level resolution.
The boundary layer itself is exactly the part already captured by the family/pair coordinates
(trailing structure) — the framework was already holding all the non-random information.

Related: [[15-borrow-chain-markov]], [[12-reload-independence]]
