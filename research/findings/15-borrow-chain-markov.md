# F15 — The ÷2-in-Ternary Borrow Chain is a (⅓, ⅔) Markov Chain

**Status: PROVEN (transition law) + MEASURED on 6,341,437 transitions from real orbits.**

## The mechanism (the problem's true home, per F2)
Halving an even number digit-by-digit in base 3 (MSB→LSB) carries a borrow state r ∈ {0,1}:
out-digit = (3r+d)÷2, new r = (3r+d) mod 2. For uniform digits:

> **P(0→1) = ⅓ (only d=1), P(1→1) = ⅔ (d=0 or 2), stationary = ½.**

## Measured on actual Collatz halving cascades
P(0→1) = 0.3513 (theory .3333), P(1→1) = 0.6458 (theory .6667).
Close but with a systematic ~2% deviation — orbit values are 3n+1 outputs, not uniform;
the residual structure is real and small. (Candidate future probe: does the deviation decay
with size, or persist? Persistent bias would be exploitable structure.)

## Conjecture E★ (the load-bearing statement, formalized)
> Along EVERY Collatz orbit, the empirical means of the index stream satisfy
> lim sup (1/m)Σ(kᵢ·log₂3 − kᵢ − wᵢ) < 0.

E★ ⟹ H decreases linearly along every orbit ⟹ no divergence. F12 shows E★ holds statistically
to four decimals; the borrow-chain Markov law is its local mechanism. Proving E★ for all orbits
(not almost all) is the program's single remaining analytic target for the divergence half.

Related: [[02-rewriting-theorem]], [[12-reload-independence]]
