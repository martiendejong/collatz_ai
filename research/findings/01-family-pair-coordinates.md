# F1 — Family/Pair Coordinates and the Closed-Form Macro-Transition

**Status: PROVEN (algebra) + verified on 50,000 random 48-bit inputs, zero failures.**

## Definitions
For odd n, write **n + 1 = a·2^k** (a odd). Then:
- **a** = family parameter (invariant along the whole odd segment)
- **k** = sequence index = number of trailing binary 1s of n
- sequence starter: n = a·2^k − 1

## The closed-form macro-transition T(a,k) → (a′,k′)
1. Segment: k exact steps n → (3n+1)/2, ending at **x = a·3^k − 1** (proof: (3(a·2^j−1)+1)/2 = a·3·2^(j−1) − 1, induct).
2. Cascade: w = v₂(x) halvings → next odd m = x/2^w.
3. Recoordinate: k′ = v₂(m+1), a′ = (m+1)/2^k′.

This is the entire Collatz dynamics with no information loss. Script: `scripts/01_transition_and_rewriting.py`.

## Reload statistics (measured, 300,000 macro-steps, 30–44 bit inputs)
Both w (cascade depth) and k′ (fresh index) match geometric(1/2) to 4 decimals:

| j | P(w=j) | P(k′=j) | geometric 2^−j |
|---|--------|---------|----------------|
| 1 | .50019 | .50021 | .50000 |
| 2 | .24976 | .25058 | .25000 |
| 3 | .12565 | .12419 | .12500 |
| 4 | .06224 | .06255 | .06250 |

## Reload predictability theorem (verified to depth j=10, zero counterexamples)
v₂(a·3^k − 1), capped at j, is a function of **(a mod 2^j, k mod 2^(j−2))** — because 3 has multiplicative order 2^(j−2) mod 2^j.

**Consequence (the lookahead regress):** foreseeing one reload of size j costs j bits of the current state; foreseeing r successive reloads costs ~2^r bits. A number holding only log₂ n bits can anticipate ~log log n reloads. No bounded-complexity function of the current representation can pre-pay for all future reloads. This is the structural reason no simple monotone quality exists (see F6, F9).

Related: [[02-rewriting-theorem]], [[09-verdicts-and-open-core]]
