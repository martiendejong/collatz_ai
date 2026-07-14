# F8 — Cycles: the circuit equation, convergents of log₂3, and a live-derived length bound

**Status: PROVEN (search ranges) + derived bound consistent with literature.**

## The circuit equation (user-conjectured, formalized, historically Steiner 1977)
A single-segment cycle ("multiply by two w times returns the ternary-converted word"):

> **a·(2^(k+w) − 3^k) = 2^w − 1**

Exhaustive search k, w < 400: unique solution **(a,k,w) = (1,1,1)** — the trivial cycle 1→4→2→1. R. Steiner (1977) proved uniqueness for all k, w via Baker's bounds; the binary↔ternary "conversion" framing is exactly his equation.

## General cycles: rational approximation rigidity
Any cycle with K odd steps, D halvings, elements > 2^71 (verification floor) satisfies
0 < D − K·log₂3 < K·log₂e/(3·2^71) → |log₂3 − D/K| < **2.04e−22**.

Continued fraction of log₂3 = [1; 1,1,2,2,3,1,5,2,23,2,2,1,1,55,…]; convergents (D/K):
2/1, 3/2, 8/5, 19/12, 65/41, 84/53, 485/306, 1054/665, 24727/15601, …, 301994/190537, 17087915/10781274, 85137581/53715833, … (Eliahou's quantization numbers reproduced exactly).

First convergent achieving error < 2.04e−22: **D/K = 103,768,467,013 / 65,470,613,321**.

> **Live-derived bound: any nontrivial cycle has K ≥ 6.55e10 odd steps, total length ≥ 1.69e11 Collatz steps.** Consistent with published bounds (Eliahou 1993 method + modern verification floor; cf. Simons–de Weger m-cycle exclusions m ≤ ~77, Hercher 2023+).

## Status of the ternary-balance intuition
"A cycle violates the ternary rule" — false as impossibility (the trivial cycle balances pen/eraser exactly; 3n−1's three cycles balance identically). True as **pricing**: the balance must approximate log₂3 to precision 1/(3·min-element), and transcendence bounds convert that into the astronomical length bound above. The rule does not forbid cycles; it makes them exponentially expensive.

## Open frontier
Multi-segment cycles with hundreds of circuits: linear-form bounds degrade with term count; the region K > 6.5e10 with m > ~77 circuits is genuinely open. This — and only this — is where a nontrivial cycle could live.

Related: [[03-pair-merge-law]], [[09-verdicts-and-open-core]]
