# collatz_ai

Autonomous AI research program on the Collatz (3x+1) conjecture — built on the
family/pair framework of Martien de Jong, executed by Jengo (Claude).

## Headline results

- **Density record:** π₁(x) ≥ x^0.9146 via the Krasikov–Lagarias inequalities at
  depths 3^13–3^20 (1.16 billion constraints at k=20, exact-integer verified), with exact-integer verified certificates (previous published
  record: 0.84, since 2003). Paper draft: `research/NOTE_DENSITY.tex`.
- **16 theorems** on the family/pair structure: Rewriting Theorem (binary
  ⟨a⟩1^k → ternary ⟨a−1⟩2^k), Pair-Merge Law, Sign Theorem + Catalan cycle
  census (complete through period 10), Shadowing Theorems, Residue-Blind
  Impossibility, exact lattice identity. See `research/NOTE.md`.
- **The Conversion Thesis:** the 3x+1 problem is the mutual illegibility of
  base 2 and base 3; binary→ternary is free, ternary→binary is the problem.
- **The eigen-visit bridge (r = 0.985):** the K–L LP-certificate's fine
  structure equals the forward orbit-visit measure (time-reversal duality).
- **Doorway anatomy:** every convergent orbit exits through an alternator gate
  (4^y−1)/3; spring gates are provably shut; 93% of all traffic passes the
  turnstiles 13 and 53.

## The current research frontier (July 2026)

- **Theorems 17-29, Lemmas 24/26/30, Propositions 23/28/31-35** in `research/NOTE.md`:
  cycle census to period 12; hop-tax exactness; zero-storage; carry characterization
  (3n+1 is the unique base-3-local map); perfect-dice factorization; CST reduced to
  rational approximation of log2(3) and verified for tau <= 4700; the master constant
  delta = log2(16/9); the drift chain (balance identity + directional low-pass +
  attenuation constant kappa_inf = 0.839).
- **Five papers** in `research/papers/`: density record (x^0.9146), tempering law,
  CST & the limma, the forgetful machine, the drift low-pass.
- **13 pre-registered predictions** for k=21 in `research/PREDICTIONS.md`.
- The gamma->1 program: every link proved or empirically pinned; remaining work is
  the variational write-up (LP duality over the K-L fixed point).

## Repository layout

| Path | Contents |
|---|---|
| `research/REPORT.html` | **The master document** — 308 research rounds, all findings |
| `research/NOTE.md` | Theorem collection with proofs |
| `research/NOTE_DENSITY.tex` / `.md` | Submission-ready density paper |
| `research/VERIFICATION.md` | 6-level anti-hallucination verification protocol |
| `research/findings/` | Per-campaign finding notes (F1–F38) |
| `research/scripts/` | All experiments (01–40), reproducible |
| `research/certificates/` | k=13/15 certificates + standalone exact verifier (k=17/19 excluded for size; regenerate via `scripts/22_kl_exact.py`) |
| `research/graph/` | Bottom-up backward graph (79,545 nodes, GraphML/GEXF) |
| `chat.txt`, `First_100_Collatz_families.csv`, … | Original source material of the family/pair framework |

## Verification

Run `python research/certificates/verify_certificates.py` — exact integer
arithmetic, no floating point in the verdict, ~2 minutes for k=13/15.

References: Krasikov & Lagarias, *Bounds for the 3x+1 problem using difference
inequalities* (arXiv math/0205002); Lagarias, *The Ultimate Challenge: The 3x+1
Problem* (AMS 2010).
