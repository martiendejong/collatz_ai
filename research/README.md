# Collatz Research Program — Family/Pair Framework

Autonomous research dossier produced by Jengo (task force session 2026-07-12), building on the
family/pair coordinate discovery in `E:\projects\collatz` (chat.txt, First_100_Collatz_families.csv).

## Structure
- `lib/collatz_lib.py` — shared machinery: (a,k) coordinates, macro-step, controls (3n−1, 5n+1), base tools
- `scripts/01…05_*.py` — reproducible experiments (each prints JSON summary)
- `findings/F1…F9` — one markdown file per finding, status-tagged (PROVEN / MEASURED / OPEN)
- `figures/*.png` — charts (funnel, streaks, needle)
- `REPORT.html` — the complete dossier in one page

## Headline results
1. **F2 Rewriting theorem** (proven): binary a·1^k → ternary (a−1)·2^k. The 111₂→222₃ mystery closed.
2. **F3 Pair-merge law** (proven + 3,002/3,002): consecutive-sequence endpoints obey x′ = 3x+2 ⇒ merge; unique to multiplier 3 ((c−1)/2 = 1 ⇔ c=3); 5n+1 misses by exactly 1 forever.
3. **F7 Needle**: never-drop survivors thin as 2^(−0.05t); survivors are trailing-ones rich (4.94 vs 2.00).
4. **F8 Cycle bound** (derived live): any nontrivial cycle ≥ 1.69e11 steps (convergents of log₂3 + 2^71 floor).
5. **F9 Open core**: measure-zero ≠ empty; Lyapunov equivalence; lookahead regress (2^r bits for r reloads).

## Controls doctrine
Every candidate argument must fail for 5n+1 (drift) AND for 3n−1 (cycle count) AND must not abolish {1,4,2}.

## Status
Conjecture: OPEN. This dossier verifies the framework, isolates the irreducible difficulty,
and ranks four concrete open leads (F9 §Open leads).
