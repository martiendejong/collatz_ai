# F22 — Period-4 Census Empty; Stopping-Time Records Reproduced Exactly

**Status: PROVEN (census) + MEASURED (records, tail).**

## Period-4 census (kᵢ, wᵢ ≤ 7; 5,764,801 shapes, sign-free)
Two integer realizations, both period-2 repeats of the −17 cycle. **Zero new.**
The integer cycle list on ℤ — {1, −5, −17} (+ boundary −1) — is now complete through period 4.

## Stopping-time records below 5×10⁶ (T-steps until value < start)
Framework hunt (restricted to n ≡ 3 mod 4 by F7's k=1 theorem):

> 27 (59), 703 (81), 10087 (105), 35655 (135), 270271 (164), 362343 (165),
> 381727 (173), 626331 (176), 1027431 (183), 1126015 (224)

This is **exactly the classical record-holder list** (Roosendaal et al.) — every famous name
recovered blind. The champions are the deepest needle-survivors: their index streams hug the
0.631 line the longest before surrendering to the drift.

## The E★ tail exponent, measured
P(stopping time ≥ t) decays at **2^(−0.080·t)** (400,000 orbits, t = 20…140).
Faster than the endpoint-only needle bound 2^(−0.050·t) from F7, as it must be: the all-prefix
constraint (never drop, not just end high) is strictly harder — the ballot-problem correction,
now measured. E★ in numbers: survival of the drift has an exponential price of 0.08 bits/step,
paid without exception in all data; the conjecture is that no orbit gets an infinite loan.

Related: [[21-period3-null]], [[07-never-drop-needle]], [[15-borrow-chain-markov]]
