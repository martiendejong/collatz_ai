# The Open Lemma, reduced to chain-flow decay

**Companion to papers/gamma_to_one.tex §Open Lemma and damping-theorem.md.
Status per component marked. 2026-07-24.**

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

## What a full proof still needs

(i) Replace measured Lemma D by a proof. Candidate route: feed domination
needs v(4m) small (Lemma B); v(4m) small needs its OWN equation to be
feed-poor AND backbone-suppressed (Lemma A prototype) — a recursive
desert-structure whose 3-adic density cost per chain level is a computable
constant < 1. The desert theorem (v₃(m+1)-cascade, Obs 319–320) is the
combinatorial engine; what must be shown is that suppression-depth requirements
accumulate additively along the chain while density decays 3-adically.
(ii) The routine-but-unwritten variance bookkeeping of the summation step
(pair-tree prefix measure, cf. damping-theorem.md Lemma 2).
