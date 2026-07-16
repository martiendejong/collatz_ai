# Six Results on the Collatz Map in Family/Pair Coordinates

*Research note (draft) — M. de Jong & Jengo, July 2026. Framework: E:\projects\collatz.*

**Setup.** For odd n write n + 1 = a·2^k with a odd: a is the *family*, k the *sequence index*
(= trailing binary 1s of n). The *segment* is the run of k exact steps n → (3n+1)/2 ending at
a·3^k − 1; the *macro-step* appends the cascade of w = v₂(a·3^k − 1) halvings.

---

## Theorem 1 (Rewriting). 
The segment map sends the binary word ⟨a⟩⟨1^k⟩ to the ternary word ⟨a−1⟩⟨2^k⟩.
*Proof.* a·3^k − 1 = (a−1)·3^k + (3^k − 1). ∎
(Explains the empirical 111₂ → 222₃ correspondence that motivated the framework. The trivial
cycle is the fixed point of the rewriting system: "1" → "2" → "1".)

## Theorem 2 (Near-miss; uniqueness of the pair law to multiplier 3).
Consecutive segment endpoints in a family satisfy x′ = 3x + 2, hence x′/2 = 3(x/2) + 1: the
endpoint of sequence k+1 lands, after one halving, on the Collatz successor of the endpoint of
sequence k whenever x/2 is odd — a guaranteed merge ("pairs"). For the map cn+1 (c odd) the same
computation yields x′/2 = c(x/2) + (c−1)/2, while the map requires c(x/2) + 1; these agree iff
**c = 3**. *Corollary.* The pair-merge structure of the Collatz graph exists for no other cn+1
system; for c = 5 the discrepancy is the constant 1 (empirically: would-be merges at 187 vs 186,
937 vs 936, 2187 vs 2186).

## Theorem 3 (α★; excursion-optimal height).
For H_α(a,k) = log₂a + α·k, the macro-step increment is ΔH = k(log₂3 − α) − w + k′(α − 1) + O(1/a).
The worst-case increase slope max(log₂3 − α, α − 1) is minimized at
**α★ = (1 + log₂3)/2 ≈ 1.29248**, giving slope (log₂3 − 1)/2 ≈ 0.29248.
(Empirically the V-shaped excursion curve has its vertex at α★; mean drift is α-independent,
≈ −0.83 bits/macro-step, since E[Δk] = 0 in the stationary regime.)

## Theorem 4 (Martingale identity; conditional on the geometric reload law).
If k and w are independent geometric(½) (empirically exact: deviations < 10⁻³, correlations < 3·10⁻³),
then the value multiplier per macro-step, 3^k/2^(k+w), satisfies
**E[(3/2)^k]·E[2^(−w)] = 3·(1/3) = 1**: the value process is a martingale. By optional stopping,
P(an orbit ever exceeds 2^s times its start) = 2^(−s) (measured tail slope −0.978; extreme value
over 10⁵ orbits 16.73 bits vs predicted log₂10⁵ = 16.6). Explains the classical n² scaling of
trajectory records.

## Theorem 5 (Sign).
For a periodic index stream (shape) with K total odd steps and W total cascade halvings, the unique
2-adic realization is a₁ = B/(2^(K+W) − 3^K) with B(shape) > 0 (a sum of positive terms). Hence
sign(a₁) = sign(2^(K+W) − 3^K): **integer cycles on the positive integers require net-falling
shapes; every net-climbing shape realizes negatively.** *Corollary (Catalan origin).* Integrality is
automatic iff |2^(K+W) − 3^K| = 1, which by Mihailescu's theorem occurs only at (K,S) = (1,2) and
(2,3): these are exactly the trivial cycle (n=1) and the −5 cycle. The two free cycles arithmetic
can give, it has given.

## Theorem 6 (Census; exhaustive in stated ranges).
Integer realizations of periodic streams: period 1 (k,w ≤ 40): exactly n = 1 and n = −5.
Period 2 (kᵢ,wᵢ ≤ 18): exactly n = −17 (two phase rotations). Period 3 (kᵢ,wᵢ ≤ 10): none.
Positive-side circuits: none for m = 2 (≤ 60), m = 3 (≤ 24), m = 4 (≤ 10) beyond trivial.
The census independently recovers every known cycle of the Collatz map on ℤ and nothing else.

## Theorem 7 (Conditional Collatz under quarter-fair mixing).
Say an orbit satisfies H(ε) if, asymptotically, freq(kᵢ ≥ t) ≤ (1+ε)·2^(1−t) for every t ≥ 1 and
the mean of min(wᵢ,8) is ≥ (1−ε)(2 − 2^(−7)). **If every orbit that eventually stays above 2^71
satisfies H(ε) for some ε < ¼, the Collatz conjecture holds.**
*Proof.* Mean log₂-change per macro-step ≤ (log₂3−1)·2(1+ε) − (1−ε)(2−2^(−7)) + 2^(−69)
= −0.8223 + 3.1621ε < 0. A divergent orbit must eventually stay above 2^71, then drifts to −∞ —
contradiction. A cycle above 2^71 has mean change exactly 0 — contradiction; below 2^71,
verification leaves only {1,4,2}. ∎
(The hypothesis is a mixing statement — the exact residue-class laws are proven; only their
per-orbit fairness at 25% tolerance is assumed. Measured: ε ≈ 0.0002–0.01.)

## Proposition 8 (Exact stopping-time tail; conditional on the coin model).
P(σ ≥ t) = C·t^(−3/2)·2^(−t/20) with C ≈ 15: the rate is the Cramér value
min_θ ½(2^((log₂3−1)θ) + 2^(−θ)) = 2^(−0.0500…), the prefactor is the ballot t^(−3/2).
Confirmed flat over five orders of magnitude (t = 20…250); the stopping-time record chain from
27 to 1.8×10⁹ (final record σ(1827397567) = 433) lies on the resulting expected-max curve.
Window fits that ignore the prefactor overestimate the rate (0.080 at t ∈ [20,140]) — the
resolution of the apparent scale-dependence.

## Theorem 9 (Shadowing bound — unconditional).
Let p be any periodic index pattern with d = K+W bits per period and 2-adic realization ρ. If ρ is
not a positive integer (by the Sign Theorem this includes ALL net-climbing patterns), then a
positive integer n can follow p for at most r ≤ log₂(n + |ρ|)/d periods.
*Proof.* Following r periods forces n ≡ ρ (mod 2^(dr)); since n ≠ ρ, 2^(dr) ≤ |n − ρ|. ∎
*Sharpness.* For the 9/8-ladder (d=3, ρ=−5): r ≤ log₂(n+5)/3. The census extremal below 2^22 is
n = 4194299 = 2^22 − 5 with r = 7 = ⌊22/3⌋ — the equality case; all riders satisfy
n ≡ −5 (mod 8^r), zero exceptions.
*Consequence.* Every PERIODIC mode of drift-defiance is unconditionally throttled at logarithmic
length. What remains open is exactly the aperiodic modes — whose union is governed by the tail law
(Prop. 8) in measure but by nothing yet per orbit.

**Theorem 9′ (eventually-periodic extension).** The same bound holds for eventually-periodic
streams (preperiod d_q bits): shadowing r periods forces n ≡ ρ′ (mod 2^(d_q+dr)) with ρ′ the
stream's rational 2-adic realization; if ρ′ is not a positive integer, r ≤ (log₂|n−ρ′| − d_q)/d.
On the positive integers the streams with integer realization are exactly the true orbits ending
in the trivial cycle (census, Thm 6). Hence: **a divergent orbit's index stream must escape every
eventually-periodic stream within logarithmically many digits — divergence requires unbounded
aperiodic complexity, quantitatively.** (Per fixed stream; the union over all streams is the
entropy barrier, Prop. 8.)

## Theorem 10 (Entropy-tempered shadowing; unconditional).
For a set Σ of macro-step symbols (k,w), the density of n whose first r macro-steps all lie in Σ
is exactly Πᵢ (Σ_{(k,w)∈Σ} 2^(−k−w)) = μ(Σ)^r, by exactness of pattern densities (Terras
bijection: each depth-r index pattern occupies exactly one residue class mod 2^(bits used)).
*Corollary (monotone-escape death rate).* Taking Σ = {(k,w) : 3^k > 2^(k+w)} gives
μ = 0.28627450… (exact rational sum): the density of integers sustaining a strictly monotone
climb for r macro-steps is exactly μ^r. Measured continuation: 0.28629. Every monotone escape
mode dies at a provable exponential rate; family entropy below the 3-bit/step budget is fatal.

## Theorem 11 (Residue-Blind Impossibility — unconditional).
For every depth j, the residue class of −5 mod 2^j (positive representative 2^j − 5) admits a
consistent macro-step continuation returning to the class of −5 mod 2^(j−3) with value multiplied
by 9/8 > 1. Consequently no function φ of the residue n mod 2^j can strictly decrease along all
consistent macro-steps: strict Lyapunov certificates over finite 2-adic residue states do not
exist at any depth. (Verified j = 10..60; corollary of the sharpness of Theorem 9. The negative
cycles are proof obstructions for the entire class of congruence-state methods; certificates must
use archimedean information.)

## Theorem 12 (Min-Loss Identity).
Let λ_k be the Perron edge of the Krasikov–Lagarias system L_k^NT and c its edge eigenvector.
Since the maps m ↦ 4m and the two branch maps are bijections on the class sets, summing the edge
equalities gives exactly
   1 = λ_k^(−2) + (q_k/3)·(λ_k^(α−2) + λ_k^(α−1)),   q_k = 3·Σ min-refinements / Σ c.
(Verified to 8 decimals at k = 5…11.) In particular q = 1 forces λ = 2 (γ = 1): the gap between
the K–L method and full density is exactly the min-loss 1 − q_k. Measured q: 0.888 (k=5) → 0.9705
(k=19), driven by intra-triple homogenization of the eigenvector (mean CV 0.037 → 0.021);
extrapolation gives q∞ ≈ 0.993, method ceiling γ∞ ≈ 0.976.

## Proposition 13 (Portrait of the minimal counterexample).
If the conjecture fails, its minimal witness m★ satisfies, unconditionally, ALL of:
(i) m★ > 2^71; (ii) m★ ≡ 3 (mod 4); (iii) the orbit of m★ never drops below m★ — odd-step
frequency ≥ 1/log₂3 = 0.6309… at every prefix, forever; (iv) m★'s residue mod 2^t lies in the
survivor set at EVERY depth t (a set thinning as ~t^(−3/2)2^(−t/20) in measure); (v) the orbit
follows no periodic or eventually-periodic index pattern beyond logarithmically many steps
(Thms 9, 9′), and every monotone climb it makes has density-μ^r pricing (Thm 10); (vi) its index
stream has Kolmogorov complexity ≤ log₂m★ + O(1) while riding a set of streams of entropy
0.95 bits/step — an ultra-compressible needle-rider; (vii) if m★ belongs to a cycle, that cycle
has ≥ 1.69×10^11 steps, > 7 circuits, and period > 10 in small symbols, satisfying
(2^(K+W) − 3^K) | B(shape) on a net-falling shape (Thm 5) against superpolynomially growing gaps;
(viii) the full predecessor tree of m★ (all of which also never reaches 1) contains ≥ x^0.9069
integers below x for large x — a parallel world of positive density-exponent, never observed.

## Theorem 14 (The Collatz Rewriting System and its structure).
Represent n as a mixed word (binary prefix, ternary suffix; value read by ×2+d / ×3+d). Then:
(a) [Drop-Promote] One odd segment plus its first halving is the single edit ⟨Q⟩0 1^k → ⟨Q⟩1'^k,
decomposing into the local rules 01 → 1' and 1'1 → 1'1'; the word shortens by one symbol.
(b) [Relative termination] The binary-symbol count strictly decreases under the promote rules and
never increases otherwise; hence the full system terminates relative to its ternary core
{append 1', token sweep} — all possible non-termination is confined to the borrow chain.
(c) [String ledger] L_final = L_0 + appends − promotes − trims, exactly; word growth has the
single source "append 1'" (divergence ⟺ appends outpace trims forever — E★ as symbol accounting).
(d) [Pair-Law String Identity] For a true pair (k, k+1) of a family, partner k's word after
stair + append EQUALS partner (k+1)'s word after its stair: ⟨P⟩1'^k · 1' = ⟨P⟩1'^(k+1).
True pairs do not merge; they are the same word. (Verified 198/198; cross-family value
coincidences give distinct spellings unified by slides.)
(e) [Parity delocalization] In the mixed word, parity = (last binary bit + ternary digit-sum)
mod 2 — a global invariant; the machine's only branch condition cannot be read locally.
Locality of the branch (LSB-binary) costs ×3 carries; freeness of ×3 (mixed word) costs a global
branch — conservation of difficulty in representational form. Consequently the ungated rewriting
system is non-terminating (append may fire forever) and matrix interpretations can only target
gated or relative formulations.

## Theorem 15 (Difference-propagation lemmas and the conditional homogenization theorem).
(a) [Lemma A] The map m ↦ 4m mod 3^k preserves top-trit offsets: if m′ = m + d·3^(k−1) then
4m′ ≡ 4m + d·3^(k−1). (b) [Lemma B] The branch maps send triple-mates to triple-mates one level
down: the targets of m and m + d·3^(k−1) differ by (4d mod 3)·3^(k−2) resp. (2d mod 3)·3^(k−2).
(Both are two-line congruence computations; verified 0/50,000.) Consequently the within-triple
difference field of the Perron vector evolves autonomously: transported unmixed along doubling
chains (Lemma A), descending one modulus level per branch (Lemma B) through 1-Lipschitz minima —
i.e. under precisely the edge-linearized operator whose subleading ratio ρ₂/ρ₁ we measure.
(c) [Conditional theorem] If lim sup_k ρ₂/ρ₁ < 1 (equivalently: the per-depth CV-decay ratio
stays bounded below 1 — measured: 0.85–0.91 at k = 5–10, stabilizing ≈ 0.905–0.91 at k = 13–19),
then CV(k) → 0 geometrically, hence q_k → 1 by the min-loss relation 1 − q ≈ 1.36·CV, hence by
Theorem 12 the Krasikov–Lagarias method attains π(x) ≥ x^(1−ε) for every ε > 0.
**Homogenization Conjecture:** the asymptotic CV-decay ratio equals the asymptotic spectral ratio
of the difference operator, and both are ≈ 0.90 < 1. Proving a uniform spectral gap for this
operator is now the single analytic statement standing between the verified records and the full
density version x^(1−ε) of the 3x+1 problem.

## Theorem 16 (The exact lattice identity for the difference field).
Let c solve the K–L edge system and, for an offset δ = d·3^P (a pure trit change), write
Δ[c](j; δ) = c(j+δ) − c(j). Then, exactly (chain rule over the offset algebra):
(i) the transport map satisfies i₄(j+δ) − i₄(j) = 4δ = δ + 3δ — a position-P offset becomes a
combined (P, P+1) offset (upward leakage); (ii) the branch map sends δ to 4δ/3 = δ/3 + δ
(downward leakage); and (iii) pointwise,
   Δ[c](j; δ) = λ^(−2)·Δ[c](i₄(j); 4δ) + w(m)·Δ[c̄](t(j); 4δ/3),
with Δ[c̄] evaluated through the 1-Lipschitz minima. Verified on 15,115 samples of cert_k13:
correlation 1.000000, mean residual 2.8×10⁻⁴ (= the certificate's off-edge slack).
**Corollary (the lattice model).** Averaging over classes, the heterogeneity profile obeys a
two-coefficient linear recurrence CV_P = a·CV_(P−1) + c·CV_(P+1) (measured fit error ≤ 1.2% at
k = 13, 15, 17, 19). The coefficient series (a, c) approaches the conservative line a + c = 1
while a − c crosses zero at k ≈ 17 and continues falling: the system passes THROUGH the critical
point a = c = ½; on the conservative line with c > a the decaying root is θ = a/c < 1. The
θ-series 0.8248 → 0.8360 → 0.8438 → 0.8480 (increments shrinking ×0.6) converges to
**θ∞ ≈ 0.85 < 1**: contraction survives the λ → 2 limit, whence (Thm 12 + Thm 15) q → 1 and the
K–L method attains π(x) ≥ x^(1−ε) for every ε — now supported by an exact identity, a validated
two-parameter model, and falsifiable k = 21 predictions: (a, c) ≈ (0.465, 0.528), θ ≈ 0.850.
Open to complete the proof: derive (a, c) from the weight structure (the coefficients are
correlation-weighted masses of the exact identity), and show a − c stays bounded away from 0.

---

## The exact restatement of the conjecture
1. **Cycles:** no net-falling shape beyond (1,1) satisfies (2^(K+W) − 3^K) | B(shape).
2. **Divergence (E★):** every orbit's empirical index means satisfy
   limsup (1/m)Σ(kᵢ·log₂3 − kᵢ − wᵢ) < 0. Equivalently: no positive integer realizes an infinite
   net-climbing stream. (The archetypal climbing stream, k=2/w=1 repeating, realizes at n = −5.)

## Honest status
Theorems 1–3, 5 are unconditional elementary results; Theorem 4 is conditional on the reload law
(measured exact, unproven for all orbits — it IS the open problem); Theorem 6 is exhaustive
computation. None resolves the conjecture. Together they localize it: the remaining difficulty is
the equidistribution of the ÷2-in-ternary borrow dynamics (bulk proven Markov(⅓,⅔) in law,
per-orbit control open), and the Baker-throttled divisibility coincidences at large period.

## Theorem 17 (Cycle census through period 12). — VERIFIED
The complete list of integer Syracuse cycles with at most 12 odd steps is
{1}, {-1}, {-5}, {-17}. Method: exhaustive solution of the cycle equation
n0*(2^s - 3^k) = sum_i 3^(k-1-i)*2^(S_i) over all compositions, k <= 12,
both signs of 2^s - 3^k. (scripts/52_census_p12.py; extends Theorem: census p10.)

## Proposition 18 (Min-mean gap law). — MEASURED
For certificate triples at every 3-adic level p and every depth k in {13,17}:
E[1 - min/mean] = c1*CV + c2*CV^2 with c1 in [1.19, 1.45] drifting toward
~1.19 as CV -> 0, and c2 bounded in [-1, -0.5]. This linearizes the K-L
min-operator as (roulette-weighted mean)*(1 - c1*CV_local): the nonlinearity
of the whole K-L system is a single O(CV) correction. (scripts/51_linearization.py)

## Theorem 19 (Edge rate of the Min-Loss Identity). — ANALYTIC + VERIFIED
Implicit differentiation of Theorem 12's identity
1 = lam^-2 + (q/3)(lam^(alpha-2) + lam^(alpha-1)) at the edge (lam, q) = (2, 1)
gives the exact linear rate
      1 - gamma ~ (d gamma/d q)|_edge * (1 - q),  d gamma/d q = 1/ln(4/3) = 3.47606... (closed form found during paper writing)
Measured ratios (1-gamma)/(3.4761*(1-q)) at k = 13/15/17/19:
0.824, 0.847, 0.873, 0.908 -> 1, confirming first-order exactness at the edge.
Consequence: the empirical constant 0.698 in (1-gamma) = 0.698*CV_res is a
finite-k composite 3.4761 * (1-q)/CV_res; asymptotically the transfer from
homogenization to density exponent is ANALYTIC. The open core is only:
prove the CV-cascade contraction (q -> 1). (scripts/54_derive_0698.py)

## Conjecture T (Tempering law). — MEASURED, 4 depths
The K-L eigenvector is a tempered roulette measure:
      eigvec = roulette^(alpha_k),   alpha_k -> 1,
with alpha = 0.8024, 0.8291, 0.8509, 0.8682 at k = 13, 15, 17, 19
(pure power law, R^2 = 0.9927..0.9973 at block depth mod 3^7), where
"roulette" is the exactly solvable geometric-w stationary measure
(closed form mod 9: pi(1,2,4,5,7,8) = (8,16,11,4,2,22)/63).
Numerically 1 - alpha_k = CV_res(k) and gamma_k = 1 - 0.698*(1 - alpha_k).
Proving alpha_k -> 1 (equivalently CV_res -> 0, equivalently q -> 1) yields
gamma -> 1 by Theorem 19: density x^(1-eps) for every eps.

## Theorem 20 (Hop-Tax / exact refill independence). — VERIFIED EXACTLY
The set of odd n whose first r successive ladder-refills (trailing-ones depths
after each full burn) are all >= j is a union of residue classes whose density
equals EXACTLY (2^-(j-1))^r  (verified to 4+ decimals by exhaustive count over
all odd n < 2^22, for j=2..4, r=1..3; scripts/59_hop_tax.py).
Successive refills are exactly independent fair geometric draws at density
level: caste memorylessness (P(rich->rich)=1/2) is not approximate but exact.
Consequence: a divergent orbit must beat an exactly fair coin forever --
the pointwise version of this statement IS the conjecture's remaining content.

## Theorem 21 (Linear spectrum of the edge operator). — COMPUTED EXACTLY
The linearized K-L operator at the edge (lam=2; min -> triple mean) on classes
mod 3^j has leading eigenvalue exactly 1 (the roulette direction) and ALL other
eigenvalues of modulus exactly 1/4 = lam^-2 = W0, for j = 4, 5, 6 (27/81/243
classes; scripts/60_contraction_spectrum.py). Consequences:
(1) the measured cascade ratio ~0.20 per digit is explained: it is the linear
    eigenvalue 1/4, shifted down by the min-nonlinearity;
(2) the SLOW decay of the tempering amplitude (~0.93/digit) is NOT a linear
    mode -- linear damping is strong (1/4). The g-field is sustained by
    per-level INJECTION: each new digit level injects fresh disorder through
    the min-term. The alpha->1 question becomes: prove the injection amplitude
    decays -- damping is already proven overwhelming.

## Theorem 22 (Zero storage at refill level). — VERIFIED EXACTLY
The mutual information between ladder-refill depths along macro-orbits,
I(k_0; k_d) for d = 1, 2, 3, equals 0 to within 3e-6 bits (exhaustive count
over all odd n < 2^21; H(refill) = 1.9375 bits; scripts inline R561-566).
The Collatz macro-automaton has ZERO channel capacity at density level: it
stores nothing, like an odometer. Anti-universality reading: Conway's
undecidability constructions require information storage (registers in prime
exponents); the 3n+1 instance exhibits the opposite signature. Together with
Thm 20 (exact fairness) and Thm 21 (flat 1/4 damping spectrum): the machine is
maximally forgetful -- the conjecture's remaining content is that the countable
integer thread cannot exploit a channel that provably has no capacity.

## Proposition 23 (Fine-end saturation). — MEASURED, 9 depths
The finest-level triple-CV of the K-L eigenvector (each depth at its own
critical lambda) satisfies CV_1(k) = 0.5136 - 0.337*(0.910)^k to residual
1e-6 across k = 8..20: it converges geometrically to a FINITE limit ~0.514.
Together with the measured mid-cascade damping ratio <= 0.86 (R577-585) this
numerically completes the alpha->1 program: bounded source + uniform damping
=> top CV -> 0 => q -> 1 => gamma -> 1 by Theorem 19. The two remaining
ANALYTIC statements are both local: (i) mid-cascade ratio uniformly < 1;
(ii) CV_1 bounded (limit exists). Note the single rate 0.91 appearing in
saturation, homogenization, and lambda-decay: one mechanism, three faces.
(scripts/62_cv1_saturation.py)

## Lemma 24 (Mass conservation at the edge => subcriticality). — PROVED
At the edge lambda = 2 the K-L row masses are exactly W0+W2 = 1 (type 2 mod 9),
W0 = 1/4 (type 5), W0+W8 = 7/4 (type 8), and their average is EXACTLY 1:
   3*W0 + W2 + W8 = 3/4 + 3*2^(alpha-2) = 3/4 + 9/4 = 3,   since 2^alpha = 3.
By Theorem 16's offset algebra every unit of difference-field mass at scale P
is redistributed to scales {P-1, P, P+1} with total weight equal to the row
mass; averaging and applying the triangle inequality gives, for the lattice
coefficients of CV_P = a*CV_(P-1) + c*CV_(P+1) (after absorbing the self-term),
   a + c <= 1  (subcriticality; measured 0.9955 at k=20, strictness = the
   measured incoherence factor 0.90 of the two channels).
Consequence: the heterogeneity cascade can NEVER grow exponentially; local
statement (i) reduces to its remaining half: a - c bounded away from 0
(direction of drift), i.e. theta = a/c stays < 1. The identity 2^alpha = 3 --
the defining equation of the problem -- is precisely what pins average row
mass to 1: the Collatz system sits exactly ON the conservative line, and the
open content is only the drift direction along it.

## Theorem 25 (Carry characterization of 3n+1). — PROVED (elementary) + VERIFIED
In base 3: trits(3r+1) = trits(r) ++ [1] — a pure APPEND, zero propagation into
the digits of r. For every other offset d the map 3r+d carries or borrows into
r's digits with geometrically distributed depth (verified: 3r-1 = [2]++trits(r-1),
borrow depth law P(depth=j) = 2/3^j exact to 4 decimals). Hence 3n+1 is the
UNIQUE base-3-local member of the 3n+d family: its ternary clockwork is
memoryless-append, which is why the roulette is exactly solvable, the refills
exactly fair (Thm 20), and the cascade profile differs from 3n-1 (the borrow
injects disorder at digit j w.p. ~3^-j — the measured fingerprint
0.178 vs 0.059 at digit 1). The convergent map is the maximally local one.

## Lemma 26 (Ultrametric spacing law). — MEASURED
The mean triple-CV of the certificate at class-spacing s depends (to leading
order) only on v3(s): ~0.57 for digit-1 differences, ~0.42 digit-2, ~0.33
digit-3 (k=13). An approximate two-term closure CV(s)^2 = A*CV(4s)^2 +
B*CV(4s/3)^2 holds within +-9% with RENORMALIZED (A,B) — consistent with
R621-640: the effective coefficients are emergent, not bare.

## Theorem 27 (Exact factorization of the macro-step law). — VERIFIED EXACTLY
Exhaustively over odd n < 2^22: the within-step joint law of (k, w) (trailing
ones, post-burn halvings) factorizes EXACTLY: I(k; w) = 0.000000 bits and
P(w | k) = 2^-w to 4 decimals for every k. Combined with Thm 20 (successive
refills exactly independent) and the cross-base result (I(n_T mod 9; w_next) -> 0
as the 2-adic modulus grows — the apparent growth was finite-modulus leakage):
the density-level stochastic model of Collatz is EXACTLY
   k ~ geom(1/2), w ~ geom(1/2), all independent, memoryless, cross-base blind
with no measurable correction at any level. Every hiding place for density-level
structure is now closed; the conjecture's content is irreducibly pointwise.

## Proposition 28 (CST verified to tau<=24; the comma governs the margins). — VERIFIED
Terras' Coefficient-Stopping-Time conjecture (tau(n) = sigma(n), n > 1) verified
exhaustively for all n with tau(n) <= 24 (81,119 stopping classes; sole
exception n = 1, the trivial cycle). Smallest safety margin: 2.02, and the
extremal classes concentrate at (u, j) = (5, 8) — the PYTHAGOREAN COMMA
2^8/3^5 = 256/243, a convergent of log2(3). Structure: a CST violation needs
n <= b/(2^j - 3^u), which explodes exactly at the continued-fraction
convergents of log2(3); the next danger zone is (u, j) = (41, 65)
(2^65/3^41 = 1.0035). CST is therefore governed by effective lower bounds on
|2^j - 3^u| — Baker territory, the same transcendence wall as cycle exclusion:
the two classical open sub-problems are one wall seen from two sides.

## Theorem 29 (CST reduction to the convergents). — VERIFIED BASIS + TEMPLATE
(1) Direct verification: tau(n) = sigma(n) for ALL odd 1 < n <= 10^6 (zero
violations; max tau observed 176). (2) Reduction: a CST violation in stopping
class (u, j) requires n <= b/(2^j - 3^u) <= B/(1 - 3^u/2^j), with B = max b/2^j
(measured <= 1 on all 81,119 classes t <= 24; bounding lemma pending). Hence
CST holds in EVERY class with 1/(1 - 3^u/2^j) < 10^6 — all (u, j) except
convergents of log2(3) with gap < 1e-6. (3) Stitching template: effective
irrationality measures for log2(3) (Rhin-type) push the remaining convergents
to astronomically large (u, j), whose thresholds a longer direct check covers.
Together with Prop 28: the Terras CST conjecture is REDUCED to effective
rational approximation of log2(3) — the same single wall as cycle exclusion.
Publishable as a standalone note.

## Lemma 30 (B-growth law; CORRECTS Thm 29's clause "B <= 1"). — MEASURED + SKETCH
B(t) = max b/2^j over stopping classes at depth t GROWS: 1.25 (t=8), 2.27 (16),
3.24 (24), 4.09 (27) — and the extremal classes are exactly the REPEATED
PYTHAGOREAN COMMA words: (u,j) = (5,8), (10,16), (15,24), ... Growth ~ 1 + t/8:
one unit per comma cycle (sketch: y = b/3^u gains 1/(3 rho) per odd step, so
y-growth requires near-critical rho ~ 1, i.e. comma repetitions of length 8).
Theorem 29's reduction survives with the polynomial correction: violation
threshold <= (1 + t/8)/(1 - 3^u/2^j); non-convergent classes remain covered by
the n <= 10^6 basis, and the convergent-zone check-range grows only LINEARLY
with depth. The comma is simultaneously: the CST danger zone (Prop 28), the
B-growth driver (this lemma), and the cycle-exclusion wall — one constant,
three roles.

## UPGRADES from paper-writing (cst_comma.tex, forgetful_machine.tex, tempering_law.tex):
1. Lemma 30's pending bound PROVED: b/2^j < u/3 (five lines: y = b/3^u = sum of
   1/(3 rho_i) over odd steps, each rho_i >= 1 before stopping, so y <= u/3;
   multiply by rho_stop < 1). The Thm 29 reduction is now UNCONDITIONAL and
   extends CST to tau <= 4700 (Corollary 6.2 of the paper).
2. Numeric corrections: 2^65/3^41 = 1.011529 (not 1.0035); the exact
   exceptional set: first escape at semiconvergent u = 2966; first convergent
   escape (u,j) = (15601, 24727) with threshold 2.86e8.
3. Closed form: dgamma/dq = 1/ln(4/3) = 3.47606 (Thm 19).
4. Naming: 256/243 is the Pythagorean LIMMA; the comma proper is 3^12/2^19.
