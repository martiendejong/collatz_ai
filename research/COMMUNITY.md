# Worldwide Collatz Knowledge Base (swept 2026-07-15)

## Verified state of the art (integrate into all our claims)
- Verification limit: ALL n < 2^71 converge (Barina, Jan 2025; J. Supercomputing 2025; 1335x GPU speedup; github.com/xbarin02/collatz)
- Cycles: no m-cycle with m <= 91 minima (Hercher 2023, arXiv:2201.00406); any cycle needs >= 1.375e11 odd terms; chain Steiner'77 -> Simons-de Weger -> Hercher
- Tao 2019: almost-all orbits attain almost-bounded values (LOG density; Syracuse random variables on Z/3^n, 3-adic Fourier decay, renewal processes)
- Krasikov-Lagarias 2003: >= x^0.84 (OUR RESULT: x^0.9069 at k=19; x^0.914 pending k=20 certification)
- Terras CST conjecture (tau(n)=t(n), coefficient stopping = stopping) — clean open subproblem, community target

## Sources
- ericr.nl/wondrous (Roosendaal): delay/glide/residue/class records; glide records complete to 2^60 (34 known); highest residue N=993 (1.253142144); delay records ~7-8/decade, ratio ~1.36
- Lagarias annotated bibliography: 363+ papers (arXiv math/0309224 + math/0608208; overview 2111.02635); umich.edu/~lagarias/3x+1.html
- ccchallenge.org: formalization project — 371 papers catalogued, 1 formalized, 4 in progress; Discord discord.gg/NJYYYMFBjA
- r/Collatz: "tuples/segments/walls" framework (mod-16 merge classes); Syracuse+predecessor-density work (u/GonzoMath); rational-cycle (3n+d) explorations
- OEIS: A006877 (delay record starters — matches our laggards exactly), A006884/A006885 (height records), A006577 (delay), A008884 (orbit of 27)
- mersenneforum + BOINC/yoyo@home: distributed verification history (87*2^60 in 2017 -> 2^71 2025)

## Community crank-failure catalogue (extends our shredder checklist)
(a) "reaches smaller number" asserted not proven; (b) reordered predecessor tree = "obvious" reachability;
(c) expected-drift 3/4 treated as proof; (d) cycles XOR divergence only; (e) mod-2^k sieve reduces density but never to 0;
(f) undecidability conflation (Conway/FRACTRAN + Kurtz-Simon Pi02 are about GENERALIZED maps, not 3n+1).
Shredder invariant: a proof must kill divergence AND all cycles for EVERY n — not density-1, not expectation, not residues.

## New facts for our program
- Our laggard records = A006877 confirmed; Roosendaal's "glide" = our stopping time sigma, "delay" = our d(n); "residue" = his Res(N)
- Hercher m<=91 supersedes our circuit exclusions (m<=7); our period-12 census is complementary (small-k exact, both signs)
- Barina's ctz+3^k-table algorithm mirrors our Drop-Promote/ladder-burn mechanics (independent rediscovery in HPC form)
- Tao's Syracuse random variables on Z/3^n = our roulette measure family (his uniform-stabilization = our homogenization; our tempering law refines: convergence is a POWER law roulette^alpha)

## Collatz as cellular automaton (attribution for our Prop 92 arc)
- Cloney, Goles & Vichniac 1987, Complex Systems 1: 3x+1 quasi-CA (base 2).
- Bruschi 2005, arXiv:nlin/0502061: two CAs for the 3x+1 map.
- Kari 2012, DLT/LNCS: base-6 local CA, Collatz & powers of 3/2, Mahler Z-numbers.
- Sterin & Woods 2020, arXiv:2007.06979: Collatz embeds base conversion (twin of our
  Conversion Thesis).

## LITERATURE VETTING OF OUR THREE NOVELTY CANDIDATES (2026-07-19)
1. GAMMA RECORD: published best is x^0.84 (Krasikov-Lagarias 2003, k=11,
   arXiv:math/0205002; earlier 0.81 Applegate-Lagarias k=9, 0.48
   Wirsching k=3, 0.43 Krasikov k=2). NO published improvement found
   since 2003 (searched 2026-07; the only recent counting preprint, Liu
   arXiv:2512.13760, claims 0.3227 - far below). OUR k=13..20 runs
   (0.8624/0.8805/0.8953/0.9069/0.9146, exact-integer certified) EXCEED
   the published record. -> genuine computational-record candidate;
   needs independent verification + writeup.
2. HILBERT/TROPICAL ANALYSIS of the K-L system: no literature found
   applying Hilbert-metric/topical-map/max-plus analysis to the
   Krasikov inequalities. Our Props 51, 61-70 (topical map, peeling
   lemma, top-band dichotomy ln4 drift, transport lemma) appear novel.
   (The "Hilbert hotel proof" preprint arXiv:2008.13643 is an
   unaccepted claimed proof - unrelated to Hilbert metric.)
3. GAMMA FORK: K-L 2002 explicitly HOPE gamma -> 1 ("one hopes x^{1-eps}
   can eventually be proved"). Our entropy-ceiling alternative
   (gamma_inf = H(1/log2 3) = 0.94996) and the k=25-30 decision test
   appear nowhere in the literature -> genuinely new open question.
MASTER SOURCES FOR ONGOING INTAKE: Lagarias' 3x+1 page
(websites.umich.edu/~lagarias/3x+1.html), his two annotated
bibliographies (arXiv math/0309224, math/0608208), ccchallenge.org
(formalization project), Tao 2019 (almost-all log-density), Hercher
2023 (no m-cycles m<=91), Barina verification, Sterin-Woods 2020.
