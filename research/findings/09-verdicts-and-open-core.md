# F9 — Verdicts Ledger and the Irreducible Open Core

**Status: synthesis of the whole program (session 2026-07-12).**

## Ledger of claims examined this session

| Claim | Verdict |
|---|---|
| Pairs merge; steps to merge determinate | TRUE — pair-merge law, 3,002/3,002 (F3) |
| Merging ⇒ everything funnels to a fixed set | INVALID — 5n+1 merges & diverges; 3n−1 merges into 3 trees |
| Stream count halves per round | TRUE per-tree; blind to tree count (F4) |
| Count monotone non-increasing | TRUE but universal — zero discriminating power |
| A branch set exists that everything falls into | TRUE, proven; sparse (density 0); cardinality open — |set|=1 IS the conjecture |
| Never-drop needs high sequence index | TRUE as lifelong prefix inequality K_j ≥ 0.631·j (F7) |
| No cycle can repeat (ternary rule) | FALSE as impossibility; TRUE as pricing → length ≥ 1.69e11 (F8) |
| Divergence needs infinite trailing 1s | FALSE — sawtooth divergence needs only finite reloads (5n+1 orbit of 7: max 12 trailing 1s at 48 bits) |
| Single-circuit cycles impossible except trivial | TRUE — Steiner 1977, reproduced (F8) |

## The one irreducible core
Every road taken this session — funnel rates, halving, index thresholds, ternary balance, height drift, streak caps — terminates at the same wall:

> **Measure zero is not empty.** All bad behavior (divergence, extra cycles) is confined to sets of vanishing density defined by infinite conspiracies of individually-legal finite events. Density arguments (Terras, Tao, and everything in this program) shrink the conspiracy sets exponentially but cannot annihilate them.

Two equivalent facts pin the difficulty:
1. **Lyapunov equivalence:** a strictly decreasing well-founded quality exists iff the conjecture is true (steps-to-1). Searching for "the key quality" IS the problem, not a shortcut.
2. **Lookahead regress (F1):** reloads are deterministic but reading r of them ahead costs ~2^r bits of state. Bounded-complexity bookkeeping provably cannot pre-pay for the future; the 3-adic/2-adic interaction scrambles faster than any representation can track.

## Open leads ranked (task force allocation)
1. **Positional alternation fuel (F6):** find a weighting of low-order 01-blocks that upper-bounds future ternlen streaks. The only lead with a plausible finite-check structure.
2. **α ≈ 1.2 height (F5):** explain the excursion minimum; test whether max-increase at α=1.2 is bounded by a constant on wider samples (if yes: revolutionary; expected: no).
3. **Cycle frontier (F8):** m-circuit exclusion past ~77 via the (a,k) chain formulation — the transition algebra may reduce term counts in the linear forms.
4. **Halving-in-ternary (F2):** the borrow-chain dynamics of ÷2 in base 3 — the problem's true home. Any equidistribution theorem here transfers directly.

## Honest bottom line
No step taken here proves the conjecture. What this program achieved: every classical partial result rediscovered and verified from the family/pair frame; two genuinely original micro-results (pair-law uniqueness to c=3 via the near-miss theorem; survivor trailing-ones enrichment 4.94 vs 2.0); and the open core isolated with unusual precision. The conjecture remains open — now with a map of exactly where it lives.
