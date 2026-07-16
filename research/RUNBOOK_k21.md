# k=21 certification runbook (lessons from five k=20 attempts)

1. CHECKPOINT every 10 power iterations (np.save) BEFORE any analysis step.
2. 60 iterations suffice (convergence by ~50; norm plateau ~1.0002).
3. Post-analysis MUST be chunked (constant memory); never stack full triples.
4. Peak RAM budget: vector (float32 N=3^20*3 = 4*3^20 bytes ~ 14GB at k=21!)
   -> k=21 needs float32 vector 13.9GB + chunk workspace: 64GB machine, or
   memmap the vector itself (slower but safe).
5. Monotone polish converges to min_ratio EXACTLY 1; margin must come from
   (a) the raw converged vector often already has margin at lambda itself
   (k=20: 1.94e-4 at 1.885) — CHECK FIRST; (b) if not, evaluate at bisected
   lambda' (note: at k=20, LOWER lambda' made the worst constraint WORSE —
   scan upward too).
6. Integer certificate: scale S=1e10, floor; margin needed only ~1e-6.
7. Exact verify: rational weight lower bounds (denominator 1e18), floor-1;
   chunked object arithmetic; ~2-4h per billion constraints.
8. gh account can be switched by other sessions: gh auth switch -u martiendejong
   before any push; big artifacts (>100MB) belong in .gitignore + Zenodo.
