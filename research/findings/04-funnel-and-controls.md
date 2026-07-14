# F4 — The Funnel: halving is real, and blind to the number of destinations

**Status: MEASURED (N = 20,000 starts, 300 steps, three systems).** Figure: `figures/e4_funnel.png`

## Stream survival (distinct current values among 20,000 tracked starts)

| step | 3n+1 | 3n−1 | 5n+1 |
|------|------|------|------|
| 10 | 8,547 | 8,546 | 16,507 |
| 50 | 2,030 | 2,111 | 13,898 |
| 100 | 579 | 558 | 13,261 |
| 200 | 71 | 51 | 13,064 |
| 300 | **3** | **25** | **13,060** |

## Readings
1. **3n+1 and 3n−1 funnel identically** — same pair law, same decay curve, step for step. But 3n+1 bottoms at 3 = |{1,2,4}| (one tree) while 3n−1 bottoms at 25 = |{1,2}∪{5,…}∪{17,…}| (three trees). **The monotone stream count converges to the total attractor size — which is exactly the unknown.** Monotonicity cannot influence its own limit.
2. **5n+1 plateaus** at ~13,060: merging exists but stalls as divergent bundles climb into sparse territory. Funnel rate ≈ drift in disguise.
3. Logical ledger: stream count non-increasing = universal (any function). Persistent strict decrease = the c=3 pair law. Neither counts trees; both are consistent with divergence (a divergent tree has infinite leaves and can halve forever while climbing).

## Proven adjacent facts
- Every trajectory merges infinitely often, or reaches 1: each odd step lands on 3n+1 ≡ 4 (mod 6), which always has the second preimage 2(3n+1). (Three-line proof.)
- Branch roots = numbers that never drop below themselves; their density is 0 (Terras 1976), log-density 0 (Tao 2019). The branch set is provably sparse; finiteness is open; |set| = 1 is the conjecture.

Related: [[03-pair-merge-law]], [[07-never-drop-needle]]
