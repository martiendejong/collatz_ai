# The Open Lemma, reduced to chain-flow decay

**Companion to papers/gamma_to_one.tex §Open Lemma and damping-theorem.md.
Status per component marked. 2026-07-24.**

> **CORRECTION (same evening, Obs 360).** The branch-period-3 identity proved
> below (D1 → D3 → D2 along the backbone) makes **Lemma A vacuous for j ≥ 2**:
> feedless backbone runs have length exactly 1, always. Horizontal
> (backbone-run) desert suppression does not exist; all suppression is
> **vertical** — through feed targets whose own subtrees are thin (the
> v₃(m+1)-cascade of Obs 319–320). Consequently the proof sketch's engine
> must be replaced by a vertical desert-stack estimate, and the numerical
> agreement of the 0.886-exponent is downgraded to suggestive until the
> vertical version is derived. The chain-flow REDUCTION itself (structure
> lemma, Lemma B, Lemma D-as-measurement, summation step) is unaffected.
> Caught by our own identity — the file is kept with this notice per the
> log's correction convention.

Setting: the K–L Perron system at the feasibility edge,
v(m) = ρ⁻¹[λ⁻² v(4m) + feed(m)], feed(m) = B·v̄(r(m)) for D1/D3, 0 for D2;
v̄ = min over the three lifts. φ(m) := feed(m)/(ρ v(m)) ∈ [0,1].
The Open Lemma (flow/L² form): intra-triple variation contracts per level by
κ ≤ κ_max < 1 in the flow-weighted L² sense.

## Where variation lives and how it moves — structure [PROVED]

1. **Backbone maps triples to triples.** The three lifts {m + t·3^{k−1}} have
   backbone targets {4m + t·4·3^{k−1}} = the triple over 4m (4·3^{k−1} ≡ 3^{k−1}
   times a unit mod 3^k). The backbone m→4m is a single N-cycle
   (ord(4 mod 3^k) = 3^{k−1}), so every triple's backbone orbit visits all
   triples.
2. **Feed targets of a triple form a triple one level down.** r₁ of the lifts:
   (4m−2)/3 + t·4·3^{k−2}, and {0, 4·3^{k−2}, 8·3^{k−2}} ≡ {0, 3^{k−2}, 2·3^{k−2}}
   mod 3^{k−1}. Same for r₃ with {0, 2·3^{k−2}, 4·3^{k−2}} ≡ triple. Hence
   intra-triple variation at level k is sourced ONLY from intra-triple
   variation at level k−1, injected through feed edges with weight φ, and
   transported along the backbone cycle with geometric weights (λ⁻²/ρ)^j.
3. **Consequence.** Unrolling the fixed point along the backbone,
   the triple-variation at m is a (λ⁻²/ρ)-geometric average, over the backbone
   orbit of m, of feed-injected variations from level k−1. Weak contraction at
   m therefore requires large φ along a RUN of consecutive backbone positions,
   and iterating in depth: **weak contraction across g levels requires a
   g-chain of feed-dominated classes** (φ > 1−ε at m, at its selected feed
   lift, at that class's selected feed lift, …).

## Lemma A (D2-run identity) [PROVED — one line]

If m, 4m, …, 4^{j−1}m are all D2 (no feed), then v(m) = (λ⁻²/ρ)^j v(4^j m)
exactly. *Proof:* iterate the defining equation; feed terms vanish. ∎
(With λ⁻²/ρ ≈ 0.29 at the edge: exponential suppression; D2-runs of length j
have 3-adic density 3^{−j}. Flow of long pure-desert stretches is thus
doubly-geometrically small — the prototype of all flow bounds below.)

## Lemma B (feed domination = backbone suppression) [PROVED — identity]

1 − φ(m) = λ⁻² v(4m) / (ρ v(m)). Hence φ(m) > 1−ε ⟺ v(4m)/v(m) < ε·ρ·λ².
A g-chain of feed domination forces the backbone side to be v-suppressed
relative to the chain at every one of its g levels: the chain lives
adjacent to a g-fold stack of suppressed backbone subtrees.

## Lemma D (chain-flow decay) [MEASURED — script 188; proof open]

F_k(g) := Perron-flow carried by classes supporting a feed-dominated chain of
length ≥ g decays geometrically in g, with ratio bounded away from 1,
uniformly in k. *(Numbers from script 188 recorded in NOTE.md Obs 358.)*

## Summation step [CONDITIONAL on Lemma D — routine]

Variance propagation across one level splits over chain length:
Var_k ≤ Σ_g F_k(g) · (amplification along a g-chain, ≤ 1 per feed edge by
φ ≤ 1) · Var-injections; with F_k(g) ≤ C θ^g, θ < 1, the sum converges
uniformly in k and yields the flow-L² attenuation κ ≤ κ_max < 1 — the Open
Lemma in the form the γ→1 chain (links 1–6 of gamma_to_one.tex) consumes.

## Proof sketch for Lemma D (the desert-pair mechanism) [SKETCH; identities PROVED]

Three exact backbone identities (machine-verified at k=9, all classes):

- **Branch period 3:** along the backbone the types cycle D1 → D3 → D2
  (m ≡ 2 → 4m ≡ 8 → 16m ≡ 5 mod 9). Every third backbone class is
  structurally feedless.
- **Adjacent feed targets are affinely locked:**
  r₃(4m) = 2·r₁(m) + 1 (mod 3^{k−1}), and the next D1 target is
  r₁(64m) = 64·r₁(m) + 42 (mod 3^{k−1}).

Now follow the chain requirement. φ(m) > 1−ε at a D1 class m demands
v(4m) < ε′·v(m), ε′ = ερλ² (Lemma B). Expanding v(4m)'s own equation, this
requires BOTH v(16m) small AND feed(4m) = B₃·v̄(r₃(4m)) small — i.e., by the
affine lock, **v̄(2·r₁(m)+1) must be desert while v̄(r₁(m)) is fertile**: the
chain forces (fertile, desert) pairs (r, 2r+1) at every level. Desert
suppression by factor δ requires desert depth j(δ) ≈ log(1/δ)/log(ρ/λ⁻²)
(Lemma A), and depth-j deserts live in specific residue classes mod 3^j
(desert theorem, v₃-cascade) of density 3^{−j}. Per chain level the density
cost is therefore

    3^{−j(ε′)} = (ε′)^{log 3 / log(ρ/λ⁻²)} ≈ (3.5·ε)^{0.886}   (at the edge),

and a g-chain costs the g-th power: **F(g) ≲ (C·ε^{0.886})^g — geometric in g
for ε < c₀**, which is Lemma D with an explicit exponent.

*Numerical check:* ε = 0.1 predicts ratio ≈ (0.35)^{0.886} = 0.395; measured
0.438 (k=12) — right order, leading-term agreement. Honest gaps: (a) the
measured k-creep of the ratio (0.438 → 0.539 by k=17) is not captured by this
k-independent skeleton — it must come from the λ(k)-dependence of ε′ = ερλ²
and of the suppression rate log(ρ/λ⁻²), both of which weaken as λ → 2; the
endpoint argument (φ ≡ 3/4 at the flat limit, so domination dies for
ε < 1/4) bounds the creep away from 1 for small ε, but the two effects must
be joined quantitatively. (b) "v(16m) small too" makes the desert requirement
recursive (a desert *stack*, not a single desert) — this only shrinks density
further, so the sketch errs on the safe side. (c) Suppression from general
desert subtrees (not just pure D2-runs) needs the recursive version of
Lemma A.

## The corrected engine: the Saturation Lemma [PROVED (sketch) — Obs 363/364]

Vertical structure of the certificate field, measured (k = 13/15) and proved:

- **Fertile tower** (m ≡ −1 mod 3^d): value RISES ≈ log₂ λ^{α−1} per rung
  (measured 0.474/0.499 vs theory 0.505/0.515 — the tower monotonicity
  v(m) ≥ (B₃/ρ)v̄(r₃(m)) is near-sharp along the −1 tower). Flow share
  decays ≈ 0.45^d: rich but rare.
- **Desert tower** (m ≡ −4 mod 3^j): penalty SATURATES at ≈ 2.5 bits,
  saturation depth 3–4. **Saturation Lemma [proved]:** for j ≥ 2 the class
  is D2 and 4m is D1 with
  r₁(4m) + 4 = (16(m+4) − 54)/3, so the inherited desert depth is
  v₃(r₁+4) = min(j, 3) − 1 (j ≠ 3): the constant **54 = 2·3³** caps the
  cascade; deep deserts (j ≥ 4) transmit depth exactly 2 regardless of j,
  and the critical stratum j = 3 splits by digit coincidence (the measured
  3–4). Since 16m is D3, a desert touches only finitely many feed
  generations: the value penalty telescopes to a bounded constant. ∎

**Consequence.** No single arithmetic structure can produce unbounded
suppression, hence none can feed a domination chain at small ε. Deep minima
of the field are multi-scale coincidences — tail events of the
log-correlated field — and the chain-flow decay (Lemma D) rests on the
field-tail estimate alone.

## What a full proof still needs

(i) **The field-tail estimate** (the analytic half, and now the ONLY half):
Gaussian-type upper bound on the lower tail of log v at each scale, with
variance growth pinned by measurement (script 191) — the machinery of
log-correlated fields / branching random walks (links 4–5 of the original
program). Everything arithmetic that could obstruct it is now proved
bounded (Saturation Lemma).
(ii) The routine-but-unwritten variance bookkeeping of the summation step
(pair-tree prefix measure, cf. damping-theorem.md Lemma 2).
(iii) Constant bookkeeping to promote the Saturation Lemma sketch to full
proof (the j = 3 stratum via measure, the finitely-many-generations
telescope with explicit constants).
