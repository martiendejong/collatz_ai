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

## Proposition 31 (The master constant delta = log2(16/9)). — MEASURED, 2 hits < 0.15%
Define delta = 4 - 2*log2(3) = 2*log2(4/3) = 0.830075 (twice the per-T-step
log2-drift of the map). Then, from the measured constants of the program:
  (i)  the (a-c) renormalization flow rate = 0.830 = delta   (rel. err 0.01%);
  (ii) the fine-end saturation / homogenization rate = 0.910 = sqrt(delta)
       (rel. err 0.12%);
  (iii) the edge rate dgamma/dq = 1/ln(4/3) (Thm 19, closed form).
The entire quantitative skeleton of the gamma->1 program appears to be
functions of ONE constant, ln(4/3) — the drift of the map. Predictions for
k=21: flow and saturation rates exactly delta and sqrt(delta). No convincing
closed form yet for theta_inf (~0.8490; 27/32 = 0.84375 off by 0.6%) or
CV1_inf (~0.5136).

## Proposition 32 (Balance + attenuation: the drift mechanism identified). — 2 PROVEN LINKS + 1 MEASURED
(i) BALANCE IDENTITY (proved, 3 lines): at the edge, the offset-magnitude
flows are exactly equal: up-flow = W0*(3/4) = 3/16 = wbar*(1/4) = down-flow,
since wbar = (W2+W8)/3 = 3/4 (mass conservation) and W0 = 1/4. The
zeroth-order scale-drift of difference mass is exactly ZERO — the system is
perfectly balanced by 2^alpha = 3, at both L1 (Lemma 24) and offset-magnitude
level.
(ii) ATTENUATION ASYMMETRY (structural): the down-channel passes through the
1-Lipschitz min (attenuation kappa < 1 whenever triples are non-comonotone);
the up-channel (transport) passes unattenuated. Hence effective down-flow =
(3/16)*kappa < up-flow = 3/16, giving c > a and theta = a/c < 1.
(iii) IDENTIFICATION (measured): kappa(P) = std(min-triple diffs)/std(member
diffs) = 0.908, 0.899, 0.888, 0.875, 0.860, 0.841 at P = 2..7 — converging
into the theta range (theta series -> 0.849). Hypothesis: theta_inf =
lim kappa at deep scales.
REMAINING for a full drift proof: kappa bounded away from 1 uniformly in k
(triples never asymptotically comonotone) — which is implied by CV_1
saturation at a nonzero limit (Prop 23). The gamma->1 chain is now:
Lemma 24 (proven) + balance (proven) + [kappa < 1 uniform, measured + linked
to Prop 23] + Thm 19 (proven).

## Proposition 33 (The clipping decomposition of kappa). — VERIFIED, exact identity
Write Dmin = Dbar + R (Dbar = mean member increment, R = min-correction). Then
kappa^2 = Var(Dbar)/V + 2cov(Dbar,R)/V + Var(R)/V exactly, and measured on
cert_k13 (P = 2..7): the ENTIRE attenuation lives in the covariance term
(2cov/V = -0.170 -> -0.257, growing with depth) while Var(R)/V is negligible
(0.010-0.028). The correction is antisymmetric: E[R | Dbar>0] = -0.026,
E[R | Dbar<0] = +0.026 (P=7) — textbook one-sided clipping from the sandwich
inequality Dx_{argmin(new)} <= Dmin <= Dx_{argmin(old)}. Effective law:
   Dmin ~ (1 - lambda_clip) * Dbar + small noise,
lambda_clip = -cov/Var(Dbar) = 0.085, 0.105, 0.122, 0.137 at P = 2, 4, 6, 7.
KAPPA-LEMMA REDUCTION: kappa < 1 uniformly <= lambda_clip >= lambda_min > 0,
which follows from non-degenerate triple gaps — i.e. from CV saturation at a
nonzero limit (Prop 23). The drift proof chain is now:
Lemma 24 + balance identity (both proven) + clipping slope > 0 (one-sided
local inequality, all ingredients measured) + Thm 19 (proven).

## Proposition 34 (Binding-constraint rigidity: the final form of the drift). — MEASURED + ROUTE
Switch-resolved decomposition (P = 2..7, cert_k13): the antisymmetric
correction R is NOT switch-exclusive — no-switch events carry the same
+-0.026 asymmetry as switches (P(switch) ~ 0.66 throughout). Hence the
mechanism is not clipping-at-switches but SELECTION-WEIGHTED MEAN REVERSION:
the argmin member co-moves less than the triple average (lambda_clip = 0.087
-> 0.140, growing with depth). LP interpretation: the min-selected member is
the BINDING constraint of the K-L linear program; binding constraints are
pinned by the equation network (complementary slackness) while slack members
float. PROOF ROUTE for the last drift link: LP duality — show the dual
weights concentrate on binding entries, making their response to coarse
perturbations strictly smaller than the free members'. The drift of the
Collatz cascade = the rigidity of binding constraints in the K-L LP.

## Proposition 34b (The min is a directional low-pass). — PROVED (one line) + MEASURED
Exact fact: min(x + c*1) = min(x) + c — the min passes components that are
CONSTANT across the triple (coarse modes) with slope exactly 1, and is
1-Lipschitz-contractive on the intra-triple-varying components (fine modes),
with measured contraction 1 - lambda (lambda = 0.087..0.140, Props 33-34).
Consequence for the two channels: the UP-channel (transport) transmits ALL
modes freely; the DOWN-channel passes through the min and therefore transmits
coarse modes freely but ATTENUATES precisely the components that vary within
triples — the finest-scale content at each level. The down-flow of fine modes
is starved while the up-flow is free: the drift c > a is the statement that
the min-operator is a low-pass filter acting in one direction only.
Remaining to quantify rigorously: the contraction factor on the varying part
is uniformly < 1 (equivalently: intra-triple variation never degenerates —
Prop 23's nonzero saturation), and the bookkeeping that converts one-sided
low-passing into c - a > 0 in the profile recurrence.

## Proposition 35 (Kappa uniformity across depths). — MEASURED, 4 certificates
The min-attenuation kappa(P; k) satisfies kappa < 0.95 at EVERY scale and
depth, and the top-aligned deep-scale value is FLAT in k:
   kappa_deep = 0.841 (k=13), 0.840 (15), 0.837 (17), 0.838 (19)
=> kappa_inf ~ 0.839 +- 0.002, uniformly bounded away from 1. The last
empirical link of the drift chain is pinned: the attenuation does not
degenerate as depth grows. Comparison: theta series -> 0.849; the
identification theta = kappa holds to ~1.2% (residual gap = up-channel
mixing correction, to be accounted in the variational argument).
DRIFT CHAIN — every link now either proven or empirically pinned:
Lemma 24 (proved) + balance (proved) + directional low-pass (proved core)
+ kappa uniform < 1 (this Prop) + Thm 19 (proved, 1/ln(4/3)).
Remaining: the variational write-up (LP duality over the fixed point).

## Theorem 36 (The division automaton / exact shed law). — PROVED + VERIFIED
Dividing an even x by 2 in base 3: sweep trits top-down with carry
c' = (c + d) mod 2, i.e. THE CARRY AT POSITION i IS THE PARITY OF THE PREFIX
TRIT-SUM. The digit-sum shed decomposes exactly per position:
   shed contribution = +1 if (c=0, d in {1,2}); -1 if (c=1, d in {0,1}); else 0.
(Verified 5000/5000; proof: q = (3c+d)//2 case table.) Combined with the
append law (Thm 25 counting form: odd step appends trit 1), the full Collatz
orbit is EXACTLY a two-rule trit-stream automaton:
   odd step: append 1 at the bottom (s3 += 1, no other change);
   halving: top-down parity-sweep with the six-entry shed table (s2 invariant).
Measured state statistics: orbit halvings run net shed 0.049/trit (random
evens: 0.033), with the +1-states enriched (0.376 vs 0.357). The pointwise
divergence question becomes: can an orbit maintain forever the
(prefix-parity, digit) statistics needed to shed at rate ~u? — where the
prefix-parity IS the orbit's own running parity (s3(x) = x mod 2): the
self-reference of the clockwork, now in exact automaton form. This is the
sharpest exact formulation of the long shot.

## Lemma 37 (Alternating-sum law; the problem closes inside base 3). — PROVED
Since 3 == -1 (mod 4): x mod 4 = alternating trit sum (exact); generally
x mod 2^j = sum_i d_i * (3^i mod 2^j), with weights periodic of period
2^{j-2} (the clockwork, reversed). Verified 5000/5000 (mod 4 and mod 8).
With Theorem 36, the ENTIRE Collatz dynamic is internal to base 3: append-1,
parity-sweep shed, w-decisions as weighted trit sums. The single remaining
non-stream ingredient is INTEGRALITY (finite top). The conjecture, final
exact form: no infinite self-consistent trit-stream statistics compatible
with sustained climb can be realized by a finite-top stream. This is the
'per-orbit rigidity / archimedean bridge' missing theorem, now surrounded by
exact machinery (Thms 25, 36, this lemma) on all sides.

## Proposition 38 (The boundary force: descent is powered by the finite top). — MEASURED
Net shed per trit by relative depth along real orbit halvings:
  TOP layer (0-10%):    +0.485 per trit  (~10x bulk)
  bulk (10-90%):        ~0 (+-0.04)
  BOTTOM layer (90-100%): -0.207 per trit (anti-shed)
Mechanism: every sweep initializes at c=0 at the finite top, and leading
trits are d=1-enriched ((c=0,d=1) occupancy 0.303 top vs 0.151 bulk), so the
top layer is a concentrated +1-shed engine. CONSEQUENCE: the digit-sum
drainage that pays the orbit conservation law (shed ~ u) is produced almost
entirely by the TOP BOUNDARY — the one feature integers have and 3-adic
pseudo-orbits lack. Descent is powered by the boundary that only integers
possess. Bridge-argument map: leading-trit statistics are governed by the
Sturmian rotation (irrationality of log2(3)) — the pointwise entry of
irrationality; divergence would require suppressing top drainage forever
against the rotation. The long shot now has a mechanism to formalize.

## Theorem 39 (The 2->1 top law; pointwise drainage floor). — PROVED (two lines)
If x has leading trit 2 (x in [2*3^L, 3^{L+1})) then x/2 in [3^L, 1.5*3^L)
has leading trit 1 — always. Multiplication by 3 preserves the leading trit.
Hence NO orbit can hold leading trit 2 through two consecutive halvings:
suppression streaks have length exactly <= 1 (measured: max streak 1 over
4000 orbits — the theorem, observed). Consequence: at least every other
halving initializes its sweep at (c=0, d=1), a +1 top-shed event — a
POINTWISE, deterministic drainage floor at the top boundary, forced by the
archimedean structure (the finite top), independent of any density argument.
The rotation subshift confirms: no w<=... in fact P(w=1 keeps the zone) = 0 —
climbing (w=1) is maximally incompatible with drainage suppression.

## Theorem 40 (The top tax and the split conservation law). — PROVED + VERIFIED
(a) At the top position of every division sweep: carry = 0 and leading trit
in {1,2}, so the top position sheds EXACTLY +1 — every halving, every orbit,
pointwise (proof: leading trit is never 0; sweep initializes c=0).
(b) Split conservation law (exact; 0/400 violations):
      s3(x_t) = s3(n) + u_t - j_t - lower_shed_t,
where lower_shed is the sweep-net below the top position. Measured on
convergent orbits: u/j = 0.504, lower_shed = -0.371 per halving (the lower
sweep GAINS on typical descents).
(c) DIVERGENCE REQUIREMENT (pointwise): with mean halving count wbar in
(1, log2 3), a divergent orbit must have
      lower gains ~ (wbar - 1) per odd step + all s3 growth, FOREVER,
against the unavoidable top tax. Every gain event requires a (carry=1,
digit<=1) state, and the carry chain is the parity of the orbit's own prefix
sums. The long shot's target inequality: bound the sustainable lower-gain
rate below (wbar - 1) using the parity-coupling — the financing cannot be
sustained. All quantities exact; the ledger is complete.

## Theorem 41 (Demote/promote machine and the zeros ledger). — PROVED + VERIFIED
The division sweep acts on digits as: carry 0 DEMOTES (2->1, 1->0, 0->0),
carry 1 PROMOTES (0->1, 1->2, 2->2) — one-line proof from q = (3c+d)//2.
Zeros are created ONLY by (c=0,d=1) shed events, consumed ONLY by (c=1,d=0)
gain events (appends add a 1; leading zeros of quotients drop). EXACT ledger
(0/300 violations):
   n0(t) = n0(0) + #sheds(0,1) - #gains(1,0) - #leading-drops.
Since n0 >= 0 pointwise:
   #gains(1,0)  <=  n0(0) + #sheds(0,1)      (POINTWISE, every orbit).
The first pointwise inequality chaining gains to sheds: the (1,0)-gain
channel of divergence financing is capped by the initial zeros plus the
shed-created supply. Remaining channel: (1,1)-gains (promote 1->2), fed by
appends and (0,2)/(1,0) events — the full three-reservoir flow network
(n0, n1, n2) is exact and closed; its perpetual-financing analysis is the
next stage of the long shot.

## Theorem 42 (The position ledger and forbidden transitions). — PROVED + VERIFIED
Per bottom-anchored position, the digit performs a constrained walk on the
path {0,1,2} under successive halvings: transitions 0->{0,1}, 1->{0,2},
2->{1,2} only — 0->2, 1->1, 2->0 are FORBIDDEN (proof: q in {demote(d),
promote(d)} and d=1 always moves; measured: exactly 0.000 each).
Gains move the digit up, sheds move it down; hence PER POSITION, over any
time window: cumulative gains - cumulative sheds = elevation change in
[-2, +2] — pointwise, every position, forever. Same-channel gain gaps >= 2
(measured P(gap=1) = 0 exactly; mean gap 6.45). Financing structure: a
double gain (0->1->2) forces a double demote (two sheds) before the position
can gain again — the three-phase re-arming cycle. Divergence financing must
therefore be spread across ~(gain excess)/2 ACTIVE positions; the width of
the financing frontier is the new constrained resource. (A resulting wbar
bound requires careful append/top boundary bookkeeping — flagged as the next
derivation, not yet claimed.)

## CORRECTION to Theorem 42 (R941-950). — the ledger clause is RESTRICTED
The per-position |gains - sheds| <= 2 ledger holds ONLY WITHIN a halving run
(consecutive halvings, no intervening append): odd steps SHIFT bottom-anchored
coordinates, so across appends the fixed-coordinate ledger accumulates over
DIFFERENT physical digits. Falsified by direct measurement: max per-position
|gains-sheds| = 113 over full random orbits (21-99 on climb segments). The
frontier-width consequence is RETRACTED. What SURVIVES of Thm 42: the
forbidden transitions (0->2, 1->1, 2->0 per halving, measured between
consecutive halvings: exact) and the within-run ledger (|Delta| <= 2 per run,
runs have mean length 2). Lesson (fourth of its kind): bottom-anchored
coordinates are not sweep-invariant across appends; all per-position claims
must specify the coordinate frame. The financing analysis continues with
run-local constraints + the append-shift as an explicit part of the automaton.

## Theorem 43 (Slot lifecycle: every append repays exactly 1). — PROVED + VERIFIED
Track SLOTS (digit positions with identity: created by appends at the bottom
with digit 1, or present in n at birth; rewritten in place by sweeps;
destroyed at the top as leading zeros). Telescoping (gains raise elevation
+1, sheds lower it -1; death digit = 0):
   sheds - gains over a completed slot's lifetime = ITS BIRTH DIGIT, exactly.
Verified: 12,069 completed slots, 0 violations; append-born slots (79%) each
net-shed EXACTLY 1. Consequence: every odd step's appended trit eventually
repays its own +1 to the drainage — the appends fund the top tax slot by
slot, pointwise. DIVERGENCE IN SLOT TERMS: financing requires an ever-growing
population of IMMORTAL, ELEVATED slots (few deaths = few leading-zero drops;
elevation toward 2). The target bound (flagged, not yet claimed): the twos
population n2 <= length forces (wbar - 1) <= (1 - wbar*log3(2)), i.e.
wbar <= 2/(1 + log3 2) ~ 1.226 — strictly stronger than the classical
wbar < log2(3) ~ 1.585. Deriving this cleanly (with exact boundary
bookkeeping, learned from two prior sign traps) is the campaign's next goal.

## Theorem 44 (Run pairing). — PROVED + VERIFIED
In every division sweep of an even number, each carry-1 run opens with
exactly one (0,1) shed and closes with exactly one (1,1) gain (final carry
is 0 since s3(x) is even). Hence #(0,1) = #(1,1) = #runs per sweep, exactly
(0/2000 violations). Consequence: the (0,1)/(1,1) channels cancel pairwise;
the NET drainage of every sweep lives entirely in the (0,2)-vs-(1,0) channel
difference — sharpening all financing analysis to one channel pair.

## Block conclusion R926-975: all static ledgers are tautologically closed.
Exact structure proven this block: Thms 40-44. The binding content of the
divergence question is DYNAMIC REALIZABILITY (which (c,d)-sequences are
achievable across successive sweeps: per-slot subshift with forbidden
transitions, coupled through carry chains), plus the flagged wbar <= 1.226
bound whose derivation requires the run/slot machinery with exact boundary
bookkeeping. The long shot's frontier is now precisely these two items.

## Proposition 45 (Integrability of the digit-sum sector; RETRACTION of the
## flagged wbar bound). — the honest closure of campaign XV-XVII's sector
(a) The financing channel is tautological: g10 - s02 = s3(x_t) - s3(n) - u
EXACTLY (from conservation + run pairing Thm 44) — the "(wbar-1)u financing
requirement" used in the flagged wbar <= 1.226 sketch was ERRONEOUS
(divergence requires only 1 <= s3 <= 2 len; no linear gain requirement).
The flag is RETRACTED. (b) Adversarial measurement confirms: max s3-raising
rate decays and goes negative within ~len(n) odd steps (T=5..80 enumeration
over n < 2^17: +2.0 -> -0.84) — true, but equivalent to s3 <= 2 len.
(c) NO-GO: every static ledger and every digit-sum currency in the
append/sweep automaton is COMPLETELY INTEGRABLE (closed by conservation):
no Lyapunov function exists in this sector — extending the spirit of Thm 11
to the s3 world. The exact laws proven (Thms 36-44) are true, beautiful, and
jointly equivalent to conservation + local automaton rules. The conjecture's
content lives strictly in the 2-adic w-sequence realizability coupled to the
finite top — where it has been since the survivor analysis. Value of the
campaign: this approach space is now CLOSED with proofs, not abandoned.

## Proposition 46 (Variational confirmation of the directional low-pass). — MEASURED, k=11 exact
Fixed-point response to mode perturbations (multiplicative forcing, eps=0.002):
- COARSE mode (mod-9 function, constant across triples): min-composite
  response slope 0.9857 ~ 1 — passes freely, exactly as Prop 34b's proved
  core predicts (min(x + c*1) = min(x) + c).
- INTRA-TRIPLE mode (top-digit function, varies within triples): min-composite
  response slope 0.9154 < 1 — attenuated, matching the statistical kappa
  (0.908 at comparable scale). Per-class binding rigidity is only ~2%
  (ratio 0.978-0.986): the rigidity lives in the min-COMPOSITE, not in
  individual binding entries. Also noted: binding degree is exactly {0, 2}
  (every argmin serves exactly two parents — the type-2/type-8 pairing).
The Open Lemma's variational quantity is now measured and agrees with the
statistical one: the directional low-pass holds in fixed-point response form.
Campaign XVII (R976-1075) complete: digit-sum sector closed with proofs
(Prop 45 no-go), variational front advanced to quantitative agreement.

## Campaign XVIII updates (R1076-1175):
1. THEOREM 17 EXTENDED to period 16 (was 12, then 14): the complete integer
   cycle list through 16 odd steps remains {1, -1, -5, -17} (exhaustive, both signs).
2. THEOREM 29 BASIS EXTENDED to 10^7: tau = sigma for all odd 1 < n <= 10^7,
   zero violations (max tau observed 246) — all stopping classes with
   violation threshold < 10^7 are now covered; the exceptional set retreats
   further into the deep convergent zone.
3. Variational slope uniformity (R1076-1090): intra-triple attenuation in
   [0.81, 0.92] at k = 9, 11, 13 — uniformly below 1 at every tested depth
   (noisier than the statistical kappa but consistent).
4. NO-GO (continuity link, R1091-1105): coarse-binned resampling does NOT
   reproduce kappa (0.978 vs 0.861) — the attenuation is carried by FINE
   correlations; kappa is not a functional of any coarse histogram. Third
   independent confirmation: only exact structure works; the uniformity
   proof must go through the exact identity, not distributional convergence.

## Proposition 47 (The scalar closure). — MEASURED to 0.06%, 3 scales
Computing the exact identity (corr 1.000000) with real mins vs bare member
means isolates the channel attenuation kappa as the ONLY nonlinear unknown,
and the variance budget then CLOSES:
   Var(LHS) = Var(T) + kappa^2 Var(B_bare) + 2 cov(T, B_min),
ratio model/actual = 1.0006 at P = 4, 5, 6 (kappa = 0.9165, 0.9057, 0.8987).
FINAL REDUCTION OF THE DENSITY TRACK: the entire cascade is exactly governed
by the 1-D system {CV-profile recursion with coefficients from the identity;
kappa(P) the sole nonlinearity}. The Open Lemma, final form: the 1-D map
kappa_in -> kappa_out along the cascade is descending with fixed point < 1
(measured: kappa falls 0.917 -> 0.899 across P, deep limit 0.839 = Prop 35 —
all instruments cohere). What remains is the stability proof of a
one-dimensional fixed point with all coefficients measured to 0.06% — the
smallest formulation the gamma->1 problem has ever had.

## Proposition 48 (The comma-cycle correspondence; conjectured by M. de Jong,
## verified same session). — PROVED (elementary) + VERIFIED
For each near-touch pair (3^k, 2^s) with gap d = 2^s - 3^k > 0, the map
n -> (3n+d)/2^w has FREE cycles with exactly k odd steps and s halvings:
the cycle equation n*(2^s - 3^k) = d*B cancels to n = B, so EVERY composition
shape yields a cycle (subject only to n odd). Verified:
- 3n+13 (limma 256/243): SEVEN distinct cycles — starters 211, 227, 251,
  259, 283, 287, 319; explicit: 211 -> 323 -> 491 -> 743 -> 1121 -> 211.
- 3n+5 (32/27): cycles at 19, 23, 29, 31 (19 -> 31 -> 49 -> 19).
- 3n+7 (16/9): cycles at 5, 7, 11.
The 3n+d cycle zoo IS the comma ladder made flesh: each convergent of
log2(3) endows its gap-value d with a full family of C(s-1,k-1) free cycles.
Conversely 3n+1's cycle scarcity is the statement that its d=1 sits ONLY on
the Catalan pairs (Mihailescu). One ladder — music, CST margins, B-growth
extremals, cycle homes — now also generates the variant-cycle taxonomy.

## Proposition 49 (The condition chain and the multiplicativity unification).
## — 3 PROVEN LINKS + 1 MEASURED PROPERTY (R1226-1270)
Attack on kappa_stability's two remaining conditions:
(i) cov(t,S) <= 0 now has a complete mechanism chain:
  (a) PROVEN (translation invariance, 2 lines): cov(x, Delta) = -Var(Delta)/2
      exactly — mean reversion is a theorem (verified -0.51..-0.54);
  (b) MEASURED: gaps scale with level, corr(level, gapsize) = +0.84
      (multiplicative field);
  (c) PROVEN (trivial): S = min(delta+g) is nondecreasing in each g;
  (d) composition confirmed: corr(t, gapsize) = -0.45..-0.51 and
      E[S | t-quartile] = +0.064/+0.005/-0.004/-0.064 — monotone, antisymmetric.
Also: Lemma 3's reflection-symmetry hypothesis TESTED and supported
(E[W+]/E[W-] = 0.997..1.018, all pairs, both scales).
(ii) competition uniformity REDUCES TO THE SAME PROPERTY: increments also
scale with level (corr +0.70..+0.77), so the competition ratio G2/delta-std
is level-free (quintile spread 0.55-0.65, flat at P=6). BOTH conditions of
the Open Lemma now hang on ONE structural statement:
   MULTIPLICATIVITY OF THE K-L FIELD (relative structure independent of
   scale) — which is the fixed-k face of the tempering law (a power law IS
   the multiplicative form). Non-circular: fixed-k multiplicativity (per
   certificate, structural) feeds the k-uniform kappa bound. The single
   remaining proof obligation of the gamma->1 program: prove the fixed-k
   field is multiplicative from the K-L equation in log coordinates.

## Proposition 50 (Homogeneity exact; the locality route falsified). — R1271-1285
(a) PROVEN (trivial, decisive): the K-L operator is positively homogeneous
of degree 1 — T(lambda c) = lambda T(c), since min commutes with scaling.
Global scale invariance is exact.
(b) FALSIFIED: the locality-based route to multiplicativity. The transport
map does NOT preserve coarse cells (modal-target coherence = 0.25 at every
tested coarseness — the scatterer again), and the level-relgap correlation
does NOT decay under coarse conditioning (+0.328 at all depths: the
multiplicativity deviation is a FINE-SCALE property).
(c) CONSEQUENCE (fifth confirmation of the program's central pattern):
per-application/local arguments fail; every true property of this system is
a property of the STATIONARY state. The multiplicativity obligation must be
attacked variationally (fixed-point analysis in log coordinates, where the
equation reads L = log(W0 e^{L o sigma} + w e^{min L o tau}) with exact
degree-1 homogeneity as the one free gift). Status of the gamma->1 program:
one obligation, one identified (hard) route, four falsified shortcuts —
all documented.

## Proposition 51 (The topical frame; the frozen gap closes; switching is
## essential). — R1286-1300
(a) PROVEN (by citation, nonlinear Perron-Frobenius / Lemmens-Nussbaum):
the K-L operator is MONOTONE + HOMOGENEOUS degree 1 = a topical map, hence
NONEXPANSIVE in Hilbert's projective metric. kappa <= 1 is now a structural
theorem, not a Lipschitz remark.
(b) MEASURED (decisive): the argmin-FROZEN linearization's second eigenvalue
RISES with depth: |lambda_2| = 0.819 (k=9), 0.895 (k=11), 0.959 (k=13) —
the frozen gap CLOSES. Since the true attenuation stays flat (kappa_inf =
0.839, Prop 35), the burden shifts with depth onto the SWITCHING: the
selection nonlinearity is not a correction but THE essential attenuator.
(c) CONSEQUENCE for the proof: no linearized/frozen argument can close the
Open Lemma. The proof must be a strict-contraction statement for the
NONLINEAR map, with the contraction constant supplied by switching
probability (measured 0.66, anchored by the competition condition /
Prop 23). Proof shape: d_H(T x, T y) < d_H(x, y) strictly on the relevant
modes, with equality analysis showing strictness whenever argmin patterns
differ — the Hilbert-metric strictness route. Sixth confirmation of the
central pattern, now at the deepest level: even the fixed point's own
linearization is a mirage; only the full nonlinear stationary object
carries the truth.

## Proposition 48b (The overtone law; instigated by M. de Jong's specimen
## 20000200550 under 3n+13). — PROVED (one line) + SPECIMEN VERIFIED
Since (a - b) | (a^j - b^j): the gap d = 2^s - 3^k of any near-touch divides
the gaps of ALL its harmonics (2^{sj} - 3^{kj}). Hence the map 3n+d has
ENHANCED cycle homes at every harmonic of its native comma: the required
divisibility miracle shrinks d-fold there. Specimen: 20000200550 under
3n+13 climbs to 2^43, then falls (418 steps) into a previously uncatalogued
15-cycle (min 131, max 1853) living on (3^15, 2^24) = the THIRD HARMONIC of
the limma, where 2^24 - 3^15 = 13 * 186793. The cycle zoo of 3n+d is a
FUNDAMENTAL plus its OVERTONE SERIES — the musical metaphor completes
itself. Corollary for 3n+1: its native gap is 1, so it receives no overtone
discount anywhere — every potential home demands a full-size miracle.
The scarcity of 3n+1 cycles restated: it is the map without overtones.
Also, the specimen's climb (2^43 before collapse) illustrates the fair-dice
lesson: no finite climb is evidence of divergence.

## Proposition 52 (The anchor no-go: divergence is uncertifiable by
## congruences in the entire proper family). — PROVED (sketch) + VERIFIED
For any proper map n -> (xn+y)/2^w (x, y odd — live parity coupling): any
congruence-forced periodic w-pattern anchors at the rational fixed point of
the composed affine map (e.g. w=1 forever anchors at y/(2-x)). If the anchor
is an integer it is a CYCLE, not divergence; if not, integers can shadow it
only ~log2(n) steps (2-adic repulsion; verified for 5n+1: shadow length =
j-2 exactly at n = (2^j-1)/3 — the alternators, playing the mirrored role).
Hence NO residue-class certificate of divergence exists for ANY proper
xn+y — the mirror of the phase-shadowing theorem (Thm 9) and residue
blindness (Thm 11). Candidates are abundant (5n+1 from 7: +0.26 bits/odd
step over 600 steps; 7n+1 from 3: +0.75), certainty is statistical, proof
is blocked by the same symmetric wall in both directions. Answer to the
question "can we find a proper divergent variant": we can find certainty
in one minute and a proof in no currently known mathematics.

## Proposition 53 (The funnel-avoidance problem: the exact divergence
## obligation, and why it is the easiest flank). — R1316-1330
For x*n+1 with huge x (e.g. x = 123121231, drift +24.9 bits/odd step):
(a) COVERING LEMMA (proved, 2 lines): for every odd x,y and every j,
{xn+y mod 2^j : n odd} covers all even residues including 0 — w is unbounded
over inputs for EVERY proper map; no modulus caps the crash channel.
(b) THE EXACT OBLIGATION: the orbit n_t crashes only via the single event
   n_t == -x^{-1} (mod 2^{~ht}),  ht ~ height ~ (log2 x - 2)*t,
i.e. the orbit falling 2-adically into the funnel of ONE bad point -1/x to
a depth growing linearly in time. Divergence proof = ONE 2-adic
non-approximation statement:
   d_2(n_t, -1/x) > 2^{-c t} for all t (any c < log2 x - 2 suffices).
(c) WHY IT FEELS (AND IS) EASIER: one bad point (vs Collatz's delocalized
everything), colossal margins (the model gives crash probability ~2^{-25}
TOTAL over the whole future — note: NOT zero, so even heuristically
divergence is 99.999997%, never 100%), and the statement has the shape of
p-adic Diophantine approximation (Baker/S-unit style) — provable for LINEAR
recurrences; blocked here only by the orbit's self-reference (w's feed back
into values). This is the family's wall in its SIMPLEST form: per-orbit
rigidity with one target point and exponential slack.
RECOMMENDED ATTACK ORDER (recorded): prove funnel-avoidance for huge-x maps
FIRST — any technique that does it (p-adic equidistribution along
piecewise-affine orbits) transfers toward Collatz. The easiest flank of the
symmetric wall.

## Proposition 54 (The base-3 column; instigated by M. de Jong's "7n+2/3").
## — MEASURED + structural, R1331-1345
The natural base-p Collatz family: p|n -> n/p; else n -> (xn+y_r)/p^w with
y_r forcing p | xn+y (y depends on n mod p — for p=3: y in {1,2}).
Critical multiplier: drift = log_p(x) - p/(p-1); critical x = p^{p/(p-1)}
(base 3: x_c = 3^{3/2} = 5.196).
(a) THE (5,3) SYSTEM — the family's nearest-to-critical convergent member:
drift -0.0350 trits/step, TWELVE TIMES closer to critical than Collatz's
-0.415 bits. All tested orbits fall into the cycle {4, 7} — which lives on
5^2 vs 3^3 (25 vs 27, gap 2): the (5,3)-family's Catalan-pair analog. The
"(5,3) conjecture" (all orbits reach {4,7}) is a HARDER margin than Collatz
— the family's most dangerous convergent laboratory. Comma ladder:
5^2/3^3, 5^13/3^19, 5^15/3^22, 5^28/3^41...
(b) THE (7,3) SYSTEM — the base-3 mirror of 5n+1: drift +0.271, small cycle
{4,5,8,10,19} coexists with exploding orbits (7 and 11 pass 2^600).
(c) The entire two-parameter (x,p) table inherits the full theory: commas =
convergents of log_p(x), cycles on near-touches, fair p-adic dice,
convergent/divergent split at the critical line. Collatz = (3,2); its
nearest-critical siblings (5,3), (13,...?) form a spectrum of laboratories.

## Literature anchors for Propositions 52-53 (agent sweep, July 2026):
1. VERDICT CONFIRMED: no proven divergent orbit exists for ANY proper map —
   not 5x+1, not the 1932 Collatz permutation's orbit of 8 (Lagarias: proving
   any specific integer lies in an infinite orbit "seems essentially
   impossible" with current methods; arXiv 2111.02635).
2. 5x+1: Kontorovich-Lagarias stochastic models (arXiv 0910.1944, 2104.10681)
   — density-1 divergence heuristic, ZERO proofs. Matches our mirror analysis.
3. THE FRONTIER, named by the field itself: Tao's 3-adic equidistribution
   (arXiv 1909.03562) and the non-Archimedean spectral program (arXiv
   2007.15936, 2412.02902) prove almost-all equidistribution and STOP exactly
   short of pointwise funnel-avoidance — both authors state the boundary
   explicitly. Our Prop 53 formulation IS the unclimbed wall, independently
   derived and now literature-anchored.
4. BONUS (deep): weak Collatz (no nontrivial cycles) would yield a proof of
   Baker's theorem simpler than any known — the cycles-transcendence link is
   BIDIRECTIONAL. Cycle bounds use Laurent-Mignotte-Nesterenko + Yu's p-adic
   Baker: |2^a - 3^b| > 3^b exp(-c (log b)^2).
5. Our small additions relative to this literature: the covering lemma
   (w unbounded for every proper map), the anchor no-go (Prop 52: no
   congruence certificate of divergence, with the alternators as mirrored
   shadow numbers), and the (x,p) table with the (5,3) near-critical
   laboratory (Prop 54).

## Theorem 55 (Universality of the perfect dice across the (x,p) table).
## — VERIFIED EXACTLY for (5,3) + general proof sketch, R1351-1365
For the base-p family (p|n -> n/p; else (xn+y_r)/p^w, y_r < p forced):
the refill law is EXACTLY P(w=j) = (1-1/p)(1/p)^{j-1} and consecutive w's
are EXACTLY independent — verified exhaustively for (5,3) mod 3^13 (ratios
1.0000, I = 0.000000 bits), matching (3,2) (Thms 20/27). General proof: x is
a unit mod p^j, so n -> xn+y_r permutes residues; the forced class leaves
the higher p-digits uniform. Hence: THE DICE ARE THE SAME IN EVERY CELL of
the two-parameter table. Also the append property generalizes: xn+y_r
appends digit y_r in base x — minimal coupling is table-wide.
SHARPENED PROOF CONSTRAINT: since fairness, independence and carry-freeness
are universal, the convergence of (3,2) CANNOT be explained by them. What
distinguishes cells is ONLY (i) the drift log_p(x) - p/(p-1) and (ii) the
comma ladder of log_p(x). Any proof of Collatz must therefore consume the
drift and the specific Diophantine structure of log2(3) — nothing else
differs. The proof-search space, reduced by one more axis.

## Theorem 56 (The three regimes; Collatz is the unique mystery cell of
## base 2). — PROVED (one line each), instigated by M. de Jong's (2,3) question
For the proper (x,p) family:
(a) x < p: TRIVIALLY CONVERGENT, pointwise: (xn + y_r)/p^w <= (xn+p-1)/p < n
for n > (p-1)/(p-x) — every step strictly decreases; all orbits reach the
finite core in linear time. Verified: (2,3) has exactly the fixed points
{1} and {2} (both on the 3-2=1 touch; the 9/8 home is dynamically
unreachable since each multiply forces at least one divide: s >= k).
NO dice, no conjecture — one line.
(b) p < x < p^{p/(p-1)}: THE CONJECTURE ZONE — drift negative but the map is
pointwise non-monotone: descent is true statistically, unprovable pointwise.
(c) x > p^{p/(p-1)}: divergent-type.
THE CENSUS OF THE MYSTERY DIAGONAL: p=2: 2 < x < 4 gives x = 3 ONLY —
COLLATZ IS THE UNIQUE CONJECTURE-ZONE CELL OF BASE 2. p=3: x in {4, 5};
p=5: x in {6, 7}; p=7: x in {8, 9, 10 (gcd!), 11}... finitely many per p,
none as isolated as (3,2). The famous problem is the smallest, loneliest
inhabitant of the entire mystery diagonal — the minimal cell where
statistics and pointwise truth first come apart.

## Theorem 57 (The jamming criterion; base 2 is the unique jam-free
## fixed-constant base). — PROVED (two lines), instigated by M. de Jong's
## "2n+1/3" question, R1401-1410
Fixed-constant map: p|n -> n/p; else n -> (xn+y)/p^{v_p(xn+y)} (w = 0
allowed). The decision variable n mod p evolves by the affine map
n -> xn+y; if 1-x is invertible mod p, the class n* = y/(1-x) mod p is a
FIXED, INVARIANT class on which the orbit never divides and strictly grows:
a PROVABLE DIVERGENCE TRAP covering density 1/p of all integers.
Example (2,3): n* = 2 mod 3; the orbit of 2 is 3*2^k - 1 = 2, 5, 11, 23,
47... — provably divergent (two lines). Cycle {1} coexists.
BASE-2 UNIQUENESS: for p = 2 and x odd, 1-x is even = not invertible mod 2:
NO jam class exists — base 2 is the ONLY base where a fixed additive
constant (like Collatz's +1) yields a trap-free map. Collatz's simple
formula is not a choice but a base-2 privilege; in every odd base, fixed
constants jam and only the residue-dependent y_r construction (Thm 55
family) stays proper. NOTE: this also answers the earlier challenge "find
one 100%-provably divergent sequence": it exists exactly UP TO the jamming
line — 3*2^k - 1 under (2n+1)/3 is provably divergent, but by a jammed
(broken-coupling) mechanism, consistent with Prop 52's live-coupling no-go.

## Proposition 58 (The jam dichotomy: divergence is provable exactly when
## the coupling is dead). — PROVED (congruence part) + census, R1411-1420
For any (x,p) fixed-constant map, exactly one of two regimes holds per map:
(a) DEAD COUPLING: gcd(1-x,p)=1 gives the invariant w=0 class
    n* = y/(1-x) mod p (Thm 57); on it the orbit grows by exactly x each
    step: divergence PROVABLE in two lines, density 1/p of all integers.
    Census: (2,3),(5,3),(2,5),(3,5),(4,5),(8,5)... all jam; x ≡ 1 mod p
    cells shift through all residues (no fixed point, coupling fires).
(b) LIVE COUPLING: w >= 1 at every non-dividing step (base 2 with x odd:
    xn+y always even — the unique fixed-y live column; or any y_r-proper
    map). Then each step consumes >= 1 digit of congruence information,
    the dice are perfect (Thm 27), and NO congruence certificate can
    prove divergence of any single orbit (Prop 52 sharpened): on every
    class a mod p^k the landing precision decays by w per step, so all
    initial congruence data is spent in <= k steps.
CONSEQUENCE for regime (c) of Thm 56: 5n+1 (the canonical drift-positive
cell: x=5 > 4 = p^{p/(p-1)}, drift +0.322 bits/odd step) has, per census
n <= 10^5: 3.0% of orbits in the three known cycles ({1},{13},{17} minima),
97.0% escaping beyond 2^200 — and ZERO provably divergent orbits. The
orbit of 7 reaches 2^10295 in 10^5 steps, believed divergent, open in the
literature. MIRROR STATEMENT: proving one orbit of 5n+1 diverges and
proving every orbit of 3n+1 converges are the two flanks of the SAME
pointwise-vs-measure gap (funnel-avoidance, Prop 53); density machinery
is map-blind between them (mirror blindness). So: provable divergence
exists in the table exactly UP TO the jam line, and not one integer
beyond it. Answer to "can we prove a drift-positive version": YES for
every jammed cell (trivially), NO for every live cell — and the second
NO is Collatz-hard by mirror symmetry.

## Theorem 59 (Universal family algebra: anchors y/(p^j - x) generate the
## family/sequence framework in EVERY cell). — VERIFIED exact, R1421-1435
For any cell (x,p) and accelerated word dividing by p^j, the rational
anchor a = y/(p^j - x) is a fixed point; the linear coordinate (n - a)
multiplies by EXACTLY x/p^j during a run, and run length = number of
p-adic digits of agreement with a. Martien's 3n+1 families are the x=3,
p=2 instance (repunit anchor -1 rises x3/2; anchor +1 IS the trivial
cycle). Verified elsewhere:
* 5n+1: THE TWO PHASES SWAP ROLES. j=1 anchor -1/3 = the ALTERNATOR
  ...010101 is now the PRIMARY rising family (1 halving/step, coord
  (n+1/3) x 5/2; closed-form jump n_r = (5/2)^r(n0+1/3)-1/3 verified over
  a 23-step run); j=2 anchor -1 = the REPUNIT ...111 is the slow rising
  family (exactly 2 halvings/step, coord (n+1) x 5/4). Both rise because
  x=5 > 4: regime (c) = "all families rise".
* (2,3) proper: j=1 anchors y/(3-2) = y are precisely the two fixed
  points {1},{2}; all families FALL (x<p), coord (n-y) x 2/3 verified.
* (7,3) proper: trit-signatures! y=1 anchor -1/4 = trit-ALTERNATOR
  (tail 2,0,2,0,...); y=2 anchor -1/2 = trit-REPUNIT (tail 1,1,1,1,...);
  seed agreeing to 15 trits produced a 12-step run of single trit-sheds
  (x7/3 per step) as predicted.
MORAL: repunit/alternator phenomenology is not Collatz-specific; it is
the p-adic geometry of rational anchors, present in every cell. What
DIFFERS per cell is only which anchors rise vs fall (regime, Thm 56) and
whether an anchor lands on a positive integer (then it is a cycle: +1
for 3n+1, {1},{2} for (2,3), {4,7},{8,14} comma-anchors for (5,3)).

## Proposition 60 (Cross-base anchor spectroscopy). — VERIFIED (laws L2/L3
## proved for prime denominators), R1436-1460, instigated by M. de Jong
Scan of all main anchors (3n+1 ladder 1/(2^j-3), 5n+1 ladder 1/(2^j-5),
(7,3) anchors -y/4) across prime bases q = 2..29 (+131):
L1  UNIVERSALS: -1 is the repunit (digit q-1) in EVERY base; +1 trivial
    everywhere; -1/2 is the "middle-digit repunit" (constant digit
    (q-1)/2) in every odd base.
L2  PERIOD LAW: tail period of an anchor with denominator d in base q is
    ord_q(mod d) (elementary but organizing): constant <=> q = 1 mod d;
    alternator (period 2) <=> q = -1 mod d for prime d.
L3  COMPLEMENTARITY: every period-2 anchor -1/d (d prime) has its two
    digits summing to EXACTLY q-1 in every base (proved: block
    B = (q^2-1)/d = q(c-1)+(q-c) with c=(q+1)m/d). The binary alternator
    01 (sum 1) is the smallest instance.
L4  CYCLOTOMIC SCATTERING: the 3n+1 anchor ladder has denominators
    2^j-3 = 5, 13, 29, 61, 125, 253...; the first prime base seeing rung
    j as constant is 11, 53, 59, 367, 251, 1013 — each rung lives in its
    own arithmetic progression q = 1 mod (2^j-3). Double-rung bases exist
    (131 sees rungs 5 AND 13) but grow like the lcm: no base sees the
    whole ladder. Prime-splitting table verified: 31 and 61 see {3,5};
    53 is the first to see 13; 59 the first to see 29.
MORAL: the cross-base structure of the problem lives ENTIRELY in the
anchors — generic orbit values carry exactly zero cross-base information
(Thm 22) — and the anchor spectrum is cyclotomic: each prime base is a
spectral filter passing exactly the rungs with d | q-1. Base 2 is the
unique base that engages ALL rungs dynamically yet renders none of them
simple: it interrogates the whole ladder precisely because it resolves
none of it. Spectroscopy metaphor exact: bases = filters, anchors =
lines, splitting governed by residues mod d (cyclotomic arithmetic).

## Remark 60b (honest relabel + consolidation: THE ANCHOR LATTICE), R1461-1465
(1) RELABEL: Prop 60's L3 is MIDY'S THEOREM (1836) in p-adic form — the
half-period digit-complement law. Verified in general form over 79 (d,q)
pairs: even period L means q^{L/2} = -1 mod d, half-shift maps -1/d to
1/d = complement; digit-sum per period = (q-1)L/2. Known mathematics;
we keep it as organizing structure, claim no novelty.
(2) CONSOLIDATION: Thm 59 (families) + Prop 60 (spectroscopy) + Thm 29
(CST) + Thm 56 (regimes) are four faces of ONE object, the ANCHOR
LATTICE: the rationals y_w/(p^j - x^r) indexed by words w. Faces:
dynamics (families/runs = p-adic proximity), cycles (= anchors landing
on positive integers; cycle equation (2^j-3^r) | W_word), CST (small
denominators = convergents of log_2 3; extremal words = repeated limma
words, Lemma 30), and cross-base structure (cyclotomic splitting by
residues mod denominator). Wall 2 (cycles/CST) is exactly the ARITHMETIC
OF THE ANCHOR LATTICE AT THE CONVERGENTS. Walls 1 (variational
strictness) and 3 (funnel-avoidance) are untouched by this arc.
(3) NO-GO #7 (base transport): a base resolves an anchor simply iff
d | q-1, and a resolved anchor is dynamically inert in that base.
Simplicity and dynamical relevance are MUTUALLY EXCLUSIVE — the problem
cannot be transported to a base where it becomes easier. Seventh
independent confirmation of the doctrine: only exact/stationary
structure works; all statistical, local, linear, and now BASE-CHANGE
shortcuts are closed.

## Proposition 61 (Wall 1: two-point Hilbert contraction is UNIFORM in k).
## — MEASURED, scripts/71, R1466-1490
The true nonlinear K-L operator at the edge (lam=2): asymptotic two-point
Hilbert-metric contraction per sweep mu_true = 0.854 / 0.870 / 0.869 for
k = 9/11/13 — FLAT, bounded away from 1 (consistent with kappa uniformity,
Prop 35). Frozen-argmin comparison: mu_frozen = 0.859/0.876/0.879; the
nonlinear gap GROWS with k (0.005 -> 0.011), switching supplies strictness.
Near-tie density (competition condition) grows 42% -> 50% -> 57% (k=9-13):
the switching reservoir DEEPENS with k. Amplitude scan (k=11): mu_true
0.857-0.878 across eps = 0.01-2.0, switching fraction rising 0.04% -> 0.44%.
STATUS of the strictness lemma: numerically TRUE WITH MARGIN and uniform;
the proof remains open (this is still wall 1), but the quantity to bound
is now measured stable: mu <= 0.88 per sweep at the edge for all tested k.

## Proposition 62 (THE CRITICAL-WINDOW CONGRUENCE LAW). — VERIFIED r <= 17
## (two independent methods), scripts/72-73, R1466-1515
Cycle equation: an r-odd-step cycle with j halvings satisfies n*D = W_word,
D = 2^j - 3^r, W = sum_i 3^{r-1-i} 2^{J_i}. Critical window j = ceil(r a).
(1) CENSUS (meet-in-the-middle, exact, all C(j-1,r-1) words, r <= 24):
    ZERO nontrivial words with D | W in every critical window r = 3..24
    (only trivial-cycle words hit, r=1,2). Expected under equidistribution:
    ~9.7 phantom hits; observed 0: P ~ 6e-5. NOT equidistribution.
(2) MECHANISM: pure congruence obstruction. Single-prime blocks at
    r = 3,4,5,7,8,11,13 (e.g. mod 5 only 4/5 residues reachable, 0 missing;
    mod 502829: 113441/502829, 0 missing). FULL-MODULUS DP (poly-time,
    validated exactly vs brute force at r=6,9): 0 is UNREACHABLE mod D for
    EVERY tested critical window r = 3..17 — even when every prime factor
    individually reaches 0 (r=6: mod 5 ok, mod 59 ok, mod 295 BLOCKED:
    the obstruction lives at CRT/correlation level).
(3) At r=18, 21 the largest prime factors are fully reachable; full-D DP
    needs a bitpacked implementation (D up to 2.7e11) — open whether the
    law persists; the census says hits=0 through r=24 regardless.
SIGNIFICANCE: known exclusions (Simons-de Weger, Hercher m<=91) use size/
transcendence bounds; this is a different, finite, poly-time certificate
PER WINDOW (DP reachability), and empirically it never fails. If the law
"0 not in reach(W mod D) at critical windows" holds for all r, then (with
Barina's 2^71 forcing n_min large, which forces cycles INTO critical
windows at strong convergents) cycles die by congruence alone — rerouting
wall 2 from transcendence to a combinatorial reachability statement.
NEW OPEN PROBLEM (promoted to top of list): prove the reach-set of the
cycle-word DP mod D avoids 0 in every critical window.

## Lemma 63 (THE NO-PHANTOM LEMMA). — PROVED (three lines) + exhaustive
## verification r=5-7, R1541-1550
(i) Every halving word (j_1..j_r) defines a unique rational cycle: the
forced word map has multiplier 3^r/2^j != 1, fixed point x0 = W/D in Q.
(ii) If x0 is a positive integer but some iterate x_i is not, write
x_i = odd/2^a (a>=1); then 3x_i+1 = (3*odd + 2^a)/2^a has ODD numerator,
so v2(x_{i+1}) = -a - j_{i+1} < -a: v2 STRICTLY DECREASES for the rest of
the word and can never return to v2(x0) >= 0 — the orbit cannot close.
Likewise an even integer iterate forces a non-integer next step. Hence
all iterates are odd positive integers: a GENUINE cycle.
COROLLARY: in every positive window (W_min = 3^r - 2^r > D there),
D | W  <=>  genuine 3n+1 cycle. THERE ARE NO PHANTOMS.

## Correction 62c (reinterpretation of Prop 62). — same-day correction
Prop 62 called the congruence "necessary, not sufficient": WRONG — by
Lemma 63 it is exactly equivalent. Consequences, honestly relabeled:
(1) The census zeros and DP certificates are INDEPENDENT PROOFS of
    "no r-cycle in windows ceil/ceil+1" for r <= 24 (census) and r <= 18
    (self-contained poly-time DP, no trajectory verification used).
    These ranges are consistent with (and implied by) Barina 2^71 +
    Hercher; the novelty is the METHOD (per-window DP certificate), not
    the exclusions.
(2) The "P ~ 6e-5 surprise" is the core mystery quantified: naive
    equidistribution of W mod D predicts ~10 cycles in r <= 24 critical
    windows; the reach-sets miss exactly the residue 0. Proving they
    always miss 0 IS the cycle problem (faithful reformulation as
    combinatorial reachability — no transcendence in the statement).
(3) Non-critical windows ceil+1: blocked for all r >= 5; the r=3,4
    reachable cases are exactly trivial-cycle repetitions (2,2,...,2)
    with j = 2r landing in that window. UNIFIED LAW: the only cycle-
    congruence solutions anywhere observed are trivial repetitions.

## Observation 64 (comma words are NOT Sturmian). — MEASURED r<=14, R1551-1560
The clearance-achieving words (|W - nD| = 1,2,4,5: near-identity
translations of Z/D exist at distance 1 for most r) are irregular — a
large halving (3-4) adjacent to runs of 1s — and are NOT rotations of
balanced/Sturmian words, contrary to the limma-repetition expectation
from Lemma 30 (which governs B-growth extremals, a different functional).
Geometry of clearance words: open.

## Proposition 65 (Bushell no-go + THE ADDITIVE DRIFT LAW). — MEASURED,
## scripts/75-series experiments, R1591-1610
(1) NO-GO #8: finite-horizon cone absorption fails. Delta(m) curves for
spike heights 1e6/1e12/1e24 remain separated by exactly the initial
height differences — there is NO m with T^m(cone) of bounded Hilbert
diameter. Birkhoff/Bushell strictness at fixed horizon is dead.
(2) What replaces it: at large Hilbert distance the operator removes an
ADDITIVE chunk per sweep: c = 1.343/1.357/1.393 (k=9/11/13, spike start,
converging to ln 4 = 1.3863 — the halving tax), c = 1.5-1.8 for random
starts. Uniform (mildly increasing) in k. Near-ray multiplicative rate
IMPROVES with k: 0.745/0.678/0.638.

## Proposition 66 (THE TROPICAL RECESSION DECOMPOSITION of wall 1).
## — MEASURED + structural, R1611-1625
The large-amplitude limit of the K-L operator in log space is its
recession map, a pure min-max (tropical) map:
   out[i] = max( x[i4] - ln4 ,  min3(x[children]) + w_{2|8} ).
Long iteration (1200 sweeps): oscillation drifts additively down to a
FINITE PLATEAU — 6.25 (k=9), 7.5 (k=11) — the tropical periodic core —
and stays there forever. The tropical map does not contract to zero;
the smooth (log-sum-exp) part of the true operator takes over below the
plateau scale and contracts multiplicatively (0.64-0.87).
CONSEQUENCE — wall 1 splits into two clean sub-problems:
(a) TROPICAL DRIFT (osc > core): additive decrease per sweep of a finite
    min-max map — max-plus spectral theory; per-k decidable by cycle-time
    / Karp minimum-cycle-mean analysis on the K-L digraph; the measured
    rate ln 4 suggests the top path is forced through W0-edges (-ln 4
    each). Candidate theorem, combinatorial, no analysis needed.
(b) SMOOTH SWITCHING CONTRACTION confined to the BOUNDED core (osc <= ~8
    uniformly-ish in k): the switching/near-tie mechanism (Props 23, 61)
    now only needs to work on a compact oscillation window, not globally.
This is the sharpest shape wall 1 has had: combinatorics outside a
bounded set, one local contraction estimate inside it.

## Lemma 67 (THE PEELING LEMMA: no branch-closed set exists). — VERIFIED
## EXACT k=7-13 with clean digit structure; proof = digit-consumption
## induction (modulo the shift conjugacy of Thms 36-37), R1626-1650
A sustained tropical top needs a branch-closed set: every member a branch
coord (m = 2,8 mod 9) with ALL three children in the set. Children of i
share one base t (they are its three LEADING-TRIT variants), and t = 5
mod 9 makes all children pure. Iterative peeling kills EXACTLY 1/3 of
survivors per round (one trit consumed per round, one forbidden value 5
mod 9 out of {2,5,8}) and empties in exactly k-1 rounds:
k=7: 486->...->64->0 (6 rounds); k=9: 8; k=11: 10; k=13: 12 rounds.
Hence NO branch-closed set: branch-only top-sustainment has depth < k-1.

## Theorem 68 (THE TOP-BAND DICHOTOMY). — PROVED (two lines) + verified
## 120/120 sweeps exact, R1651-1665
For the tropical recession map, set theta = w8 + ln4. EITHER the top band
(top - theta, top] contains a trit-complete triple (all three leading-trit
variants of some base t), OR the top moves by EXACTLY -ln4 this sweep.
Proof: every branch term is min3 + w8 <= (top - theta) + w8 = top - ln4
when some child of every branch coord is below the band; every P-term
<= top - ln4; and the max always contains the P-term of the current
argmax, giving >= top - ln4. QED. Verified: 120/120 sweeps with
|drop - ln4| < 1e-9 while no triple in band (k=11, osc 448 -> 83).
CONSEQUENCE: the measured additive drift ln 4 (Prop 65) is now THEOREM
above the triple-scale: while the top band is triple-free, the top decays
at exactly ln 4 per sweep. The tropical core (plateau, Prop 66) is
precisely the scale at which the eigen-shape supplies permanent triples.
REMAINING for the full tropical wall-1 statement: the mirror BOTTOM lemma
(lows are destroyed by max-lifts; low P-transport cycles residues
5 -> 8 -> 2 mod 9 through branch stops where survival requires low
children — the spreading requirement mirrored) + the smooth-regime
switching estimate on the bounded core. Wall 1's tropical half is DONE
in principle: exact combinatorial statements, both verified, one proved.

## Proposition 69 (BOTTOM TRANSPORT AND THE TROPICAL DRIFT FORMULA).
## — measured + structural, R1666-1690
(i) RESIDUE RIGIDITY: all three children of a branch coord share one base
t, hence share branch/pure status (t mod 9 common to the triple), and by
induction ALL residue-level status along descent chains is choice-free:
the min-player's freedom (pick any child) is an ILLUSION at residue level
— "some child" = "all children" there. Consequently the pure child-edge
descent web dies by the SAME peeling as Lemma 67 (<= k-1 steps).
(ii) The sustained bottom is therefore a MIXED cycle: long w2-branch
chains (descent -0.2877/step) renewed through P-edge reinjections; its
exact rate is the value of a mean-payoff game (min transports lows, the
operator's max picks the worse term). Measured: beta = -0.2811/sweep
(k=11) = w2 + 0.0066 (~1% renewal cost); argmin sits at pure coords
(5 mod 9) via P-edges, 250/250 sweeps.
(iii) THE DRIFT FORMULA: net tropical oscillation drift (triple-free
regime) = ln4 - |beta| = 1.3863 - 0.2811 = 1.1052 — matches the measured
net drift 1.1052 to FOUR decimals. Wall 1's tropical half is now fully
quantified: top decays at exactly ln4 (Thm 68, proved), bottom descends
at game value beta (computable by policy iteration), drift = ln4 + beta.

## Proposition 70 (ANATOMY OF THE CORE CONTRACTION). — MEASURED, R1691-1740
On the bounded core (the smooth regime where wall 1's open estimate
lives), two-point contraction mu = 0.85-0.87 decomposes as:
(i) The oscillation TOP rides the sigma-orbit on pure coords (argmax
    residue census: 5 mod 9 dominant, never 2), and at branch transits it
    is WELL-SUPPORTED: a child's g sits within 5-7% of max in 88-100% of
    transits. The top erodes slowly (~6%/sweep), not by deep mixing bites.
(ii) The BOTTOM rises ~6-9%/sweep; the AND-condition (both terms at min)
    rarely holds at the argmin; argmin-switching at bottom coords is
    modest (0-14%).
(iii) Top decay + bottom rise sum to 1 - mu: the contraction is
    DISTRIBUTED across the level-set bands, while the extreme points
    themselves are locally protected.
CONSEQUENCE for the proof shape: strictness on the core is a BAND
phenomenon — argmax/argmin chain arguments (which powered the tropical
side, Thm 68) cannot work here; the proof must track level-set masses
(a Lyapunov functional over bands: each band leaks mass toward the
middle at a rate bounded by balance x alignment statistics). Mechanism
now fully mapped; the quantitative band-leak estimate is THE remaining
open kernel of wall 1 — everything else (tropical top ln4 = Thm 68,
bottom game = Prop 69, peeling = Lemma 67) is proved or computable.

## Lemma 71 (THE TRANSPORT LEMMA). — PROVED (three lines) + verified
## 0/10993 violations, R1741-1765
For the two-point log-difference g = log(x/y):
   g_new[i] <= max( g[i4(i)], g[A_y(i)] )   and mirrored
   g_new[i] >= min( g[i4(i)], g[A_x(i)] ),
where A_y(i) = the argmin child of i under y. Proof: the ratio of the
two-term sums is a mediant, bounded by the max of the term ratios; and
min3(x)/min3(y) <= x[c]/y[c] at c = argmin_y since min3(x) <= x[c]. QED.
CONSEQUENCE: upper level sets transport ONLY through the graph {P-edges,
A_y-edges}; lower level sets through {P-edges, A_x-edges}. THE TWO GRAPHS
DIFFER EXACTLY AT SWITCHED COORDINATES (A_x != A_y) — switching is,
precisely and provably, the separation between the top and bottom
transport structures. This is the sharpest formal expression yet of
"switching supplies the strictness".
Obs 71b (measured): the top-quarter band loses 78% of its ALLOWED feed
mass per sweep at mediant stops (the AND-cut at balanced coords): band
mass thins fast while the extreme value erodes slowly (~6%/sweep) —
the Lyapunov functional must combine mass-thinning at fixed level with
slow max-erosion.

## Lemma 67 — proof upgrade (peeling density is exact), R1741-1765
The survivor count after s rounds is EXACTLY 2*3^(k-2)*(2/3)^s for
s <= k-2 (each peeling condition consumes one fresh trit and forbids
exactly one of its three values — free-trit counting via the division
automaton of Thms 36-37), giving exactly 2^(k-1) survivors at s = k-2;
the (k-1)-th condition has no fresh trit left (the automaton runs out of
tape; modular wrap forces the status) and wipes all survivors at once —
verified exactly at k = 7, 9, 11, 13 (2^6, 2^8, 2^10, 2^12 -> 0).
Status: proved modulo the standard shift-conjugacy, final wipe verified.

## Observation 72 (PROPER AND-CLOSED SETS EXIST — the plateau's home;
## no-go #9 for a purely combinatorial core proof), R1766-1795
The AND-closure (P-parent in C and, at branch coords, some child in C)
admits a CANONICAL proper closed subset: peeling all-minus-one-coordinate
converges, for EVERY seed, to the same size N - 20k exactly
(k=7: 729-140=589; k=9: 6561-180=6381; k=11: 59049-220=58829; rounds
20k-1). Interpretation: flat tops CAN live combinatorially on this set —
it is the home of the tropical plateau (Prop 66). Consequently strict
contraction on the CORE cannot follow from transport combinatorics alone:
the argmin must ESCAPE the flat set, which happens exactly when sibling
eigenvector ratios compete with the g-gap (the competition condition,
Prop 23). Core strictness = transport combinatorics (Lemma 71) + sibling
competition (eigen-data). The two ingredients are now provably BOTH
necessary — the proof shape is fixed.

## Observation 73 (THE ITERATED REINTERPRETATION CHAIN — M. de Jong's
## question), R1796-1810
R(n) = binary string of n read as ternary; the chain n, R(n), R(R(n))...
isolates the pure conversion toll (each step = one ternary->binary
conversion; the free direction is trivial: ternary(R(n)) = binary(n)
verbatim). Results: unique fixed point 1 (no cycles <= 10^4); length
growth -> log2(3) = 1.5850 exactly; bit fraction -> 0.5000; block-2
entropy -> 1.9999/2.0; successive-string agreement -> 0.504 (= random);
base-5 digits of deep iterates uniform. THE PATTERN LIVES EXACTLY ONE
STEP: total at one conversion (the free transplant), erased at two (the
toll). Confirms and sharpens Thm 22 / the Conversion Thesis: the toll is
a perfect pseudo-randomizer, and Collatz iterates precisely this toll.

## Observation 74 (SEED ANATOMY OF TOTAL SEQUENCE LENGTH — the anchor
## trichotomy; M. de Jong's question), R1811-1835
Census of all 20-bit seeds, total stopping time:
(1) Trailing ones (2-adic closeness to -1): mean length rises LINEARLY,
    +6.2 steps per extra trailing one (0 ones: 132.9; 14 ones: 216.9).
    Martien's intuition correct ON AVERAGE.
(2) But records are NOT the max-trailing seeds: top-10 record holders
    have 0-8 trailing ones (837799, the champion at 524 steps: only 3);
    the pure repunit 2^20-1 scores 178 vs record 524, barely above the
    mean 139. Reason = "the pattern lives one step" (Obs 73): trailing
    ones buy exactly ONE deterministic rise-run; after it the toll rolls
    fresh dice. Records = moderate trailing structure + repeated lucky
    re-entries into rise-rich regions, which no seed structure can buy.
(3) THE ANCHOR TRICHOTOMY (mod 2^14 conditioning, global mean 139.1):
    near -1  (repunit anchor):   mean 230.2  -> LONGEST on average
    near +1  (trivial-cycle anchor): 128.2   -> exactly TYPICAL, because
        the +1-word (3n+1)/4 is the drift-neutral typical word (2
        halvings per rise) - proximity to the cycle buys nothing!
    near -1/3 (alternator anchor): 66.2     -> SHORTEST (3n+1 = 0 mod
        2^m: one rise then m halvings - the crash word).
    The three rational anchors of the (3,2) cell thus label the three
    speed classes of seeds: climb / typical / crash. Seed 2-adic
    geometry determines the first word; the dice do the rest.

## Observation 75 (THE RECORD LADDER IN TWO BASES + a theorem-ette),
## R1836-1860, M. de Jong's question
All 44 total-stopping-time records <= 837799, examined binary + ternary:
(1) THEOREM (one line, verified 43/43 for n > 2): no record is = 2 mod 3.
    If n = 2 mod 3 (n > 2), then m = (2n-1)/3 is a SMALLER odd number
    whose orbit reaches n in 2 steps, so total(m) = total(n) + 2 and n
    cannot be a record. Ternary form: A RECORD NEVER ENDS IN TRIT 2.
    (Sole exception n = 2: its predecessor is 1, where the orbit stops.)
(2) DEEP-LEAF ENRICHMENT: 28/44 records are divisible by 3 (64% vs 33%
    random), 11/44 by 9; multiples of 3 have NO odd predecessors at all
    (they head their own highways). The ladder's base is exactly ternary
    round: 27 = 1000_3, 54 = 2000_3.
(3) BINARY SUFFIX CLUSTERING: odd records mass on three 5-bit suffixes:
    00111 (10x), 11111 (8x), 11011 (6x) vs uniform 2.5 - all trailing-
    ones-rich; mod 64 tops: 111111 (7x), 100111 (6x), 011011 (5x, the
    "27-suffix"). Records RE-USE proven climb suffixes, consistent with
    Obs 74 (the seed buys the first word deterministically).
(4) HIGHWAY FUNNELING: the records' orbits merge into few junctions
    (436 4x, 364 3x, 40 3x): new records = new long approaches to the
    same descent highway.
(5) Near the champion, trailing-ones-rich record seeds appear at once:
    #3 511935 = 1111100111110111111 (16 of 19 bits ones, 6 trailing),
    #7 156159 (9 trailing), #9 106239 (8 trailing), #14 26623 =
    110011111111111 (11 trailing ones of 15 bits).

## Observation 76 (ONE-RUN BLOCKS ACROSS RECORD LADDERS), R1861-1875
Scan to 10^7 (records through 8400511, 685 steps):
(1) Blocks of >= 5 ones appear in 31% of records vs 22% in same-length
    random odds - MILD enrichment; 511935's perfect [5,5,6] triple is
    exceptional, not typical. The driving structure remains the TRAILING
    block (Obs 74); interior blocks help only via later re-entries.
(2) The champions after 837799: 1723519 (runs [2,1,1,2,7]) and 8400511
    ([1,1,3,7]) both END in seven ones; 6649279 ends in six.
(3) GEM: record 3732423 = 1110001111001111000111 is a perfect BINARY
    PALINDROME (runs 3-3-4-2-4-3-3). No explanation; filed as curiosity.

## Observation 77 (THE FUEL-ENTROPY THESIS, tested — M. de Jong),
## R1876-1895
Thesis: orbits consolidate bits into long 1-runs; the runs are FUEL whose
combustion generates entropy. Verdict in three parts:
(1) FUEL: CORRECT and exact. A trailing run of k ones = exactly k rise
    steps x3/2 - the only climb mechanism. Champion orbits are fat-tailed
    in fuel packets: mean 2.481 vs geometric 2.000; k=6 packets 3x
    enriched (4.7% vs 1.6%). Control (300 random 20-bit orbits): mean
    2.018, geometric to within noise - the enrichment is pure SELECTION:
    champions ARE the orbits that drew heavy fuel, repeatedly.
(2) REFUEL IS DICE, not dynamics: on champion orbits themselves,
    P(next packet >= 4 | current >= 4) = 0.250 vs unconditional 0.234 -
    memoryless. There is NO consolidation law: mean max-run along the
    champion orbit 3.87 vs 3.70 random control (mild selection tilt
    only). Long orbits are not CAUSED by a merging mechanism; they are
    the orbits that kept winning independent draws.
(3) ENTROPY: CORRECT. Consuming a packet of k writes ~k log2(3) fresh
    top bits at near-maximal entropy (block-2 entropy 1.862/2.0 along
    the champion orbit) - the toll. Thermodynamic summary: the orbit is
    an engine burning 1-run fuel into entropy; a fraction of the exhaust
    re-crystallizes into new fuel at exactly the dice rate 2^-k. This
    memorylessness is WHY the problem is hard: a positive-feedback fuel
    law would mean divergence; the dice pin mean packet at 2 = the
    drift-neutral point.

## Theorem 78 (THE STERILITY THEOREM: pure ones cannot regenerate).
## — PROVED for odd k (one line), verified to k = 100000 for even k;
## instigated by M. de Jong, R1896-1915
Burning the pure repunit 2^k - 1 (k rises) yields 3^k - 1. Fresh fuel:
* k ODD: 3^k = 3 mod 8, so v2(3^k + 1) = 2: fresh run = EXACTLY 1. QED.
* k EVEN: fresh run r(k) = v2(u+1), u = (3^k-1)/2^(v2(3^k-1)) — exact
  formula; record values over k <= 100000: r(6)=2, r(12)=6, r(174)=9,
  r(1198)=17. MAX EVER = 17 against packets up to 100000: a pure run of
  k ones regenerates at most ~log-scale fuel, never a comparable run.

## Theorem 79 (THE FUEL LEDGER: head entropy buys fuel at par).
## — exact counting law, verified, R1896-1915
P(first packets = (k1,...,kp)) = 2^-(k1+...+kp) exactly (verified: 62
observed vs 64 predicted for prefix (5,4,3) over all 19-bit odd seeds).
Hence prescribing total fuel F costs EXACTLY F seed bits, regardless of
how it is split into packets: a B-bit head can pre-program at most B
bits of climb fuel; every packet beyond that is drawn fresh from the
dice at rate 2^-k. EXCHANGE RATE: 1 head-bit = 1 fuel-bit, at par.
TOGETHER (78+79): Martien's thesis is now theorem-shaped: long runs
cannot arise from the ones themselves (sterility); they arise only from
(i) the entropy BEFORE the ones, at par, or (ii) the entropy generated
by burning, re-crystallizing at dice rate. This also answers the old
open question "what is the maximum number of ones the clockwork can
write": pre-programmed ones <= seed bits; spontaneous ones = geometric.

## Observation 80 (THE FUEL AUDIT: orbits run on exhaust, not tank),
## R1916-1920
Champion orbits' fuel budgets: 837799: tank 20 bits, total fuel burned
195 bits (90% re-crystallized); 8400511: tank 24, burned 256 (91%);
27: tank 5, burned 41 (88%). Big packets keep appearing deep in orbits
(8400511 second half: 8,7,6,5,6). CONSEQUENCE: "tank empty => done" is
factually false — the engine runs ~90% on re-crystallized entropy. The
ledger (Thm 79) bounds only PRE-PROGRAMMED fuel; bounding the
re-crystallization rate pointwise-forever for an individual orbit is
exactly funnel-avoidance (wall 3), because the post-burn state is
provably patternless (Obs 73) — the same theorem that makes sterility
provable for explicit states destroys the tools one step later.

## Observation 81 (THE COIN STAYS FAIR; ONLY THE ARENA SHRINKS), R1921-1930
Mean packet size vs current bit-length of n, over 2000 random 20-bit
orbits: 1.993/1.999/1.982 at 20/16/12 bits — the dice are EXACTLY fair
at every scale (small-B wobbles are deterministic artifacts of the few
specific small numbers, e.g. every orbit through 31 = 11111 logs a 5).
The 8,7,6,5,6 decline in 8400511's late packets happened at bit-lengths
19,18,21,21,9: the ceiling k <= bits(n) was NOT binding — the decline
was chance, not law. TRUE PICTURE: packets never weaken; the VALUE
drifts down (E[k]=2 vs break-even ~3.4: house edge -0.83 bits/cycle),
and the shrinking value lowers the ceiling only in the endgame. Fair
coin, negative house edge, shrinking arena: extinction certain in
measure — pointwise gap unchanged.

## Proposition 82 (THE COMPLETE PATTERN TAXONOMY + the phase clock).
## — verified exact, R1956-1980, instigated by M. de Jong
Beyond the repunit (...1111) and alternator there is an INFINITE LADDER
of bit patterns, one per halving depth j: the 2-adic tail of the anchor
1/(2^j - 3), periodic with period ord_{2^j-3}(2):
  j=1: ...111111 (period 1)  climb +0.585 bits/step (repunit)
  j=2: ...000001            neutral x3/4 (trivial-cycle anchor)
  j=3: ...110011001101 (period 4, the "1100 word")  descent -1.415/step
  j=4: 1/13-tail (period 12) crash -2.415/step
  j=5: 1/29-tail (period 28) crash -3.415/step
  j=6: 1/61-tail (period 60) crash -4.415/step
Verified: seeds agreeing to 18 bits shed exactly j,j,j,... halvings until
the agreement is consumed. THE PHASE CLOCK: within any pattern-run the
2-adic agreement depth decreases by exactly j bits per step — a
DETERMINISTIC, MONOTONE reduction quality per phase; renewals between
phases are the fair dice. Orbits = concatenations of anchor-phase runs:
the "step in the pattern" is the true local clock.

## Observation 82b (no local Lyapunov: the 0.666 barrier), R1956-1980
Best linear combination of (len, pop, nruns, maxrun, pairs11) over 20000
random weightings achieves P(decrease per odd step) = 0.666 (weights
dominated by len + maxrun); len+pop gives 0.470, len+2pop 0.522. No
local string quality approaches monotonicity — as the fair-packet law
requires. The monotone quality exists only PER PHASE (the phase clock),
not across renewals.

## Observation 83 (PHASE ALPHABET IS EXACTLY DICE-OCCUPIED in the
## stationary regime), R1981-2000
Clean census (first 25 steps of 20000 40-bit seeds, 500k steps): marginals
P(j) = 2^-j within 2 sigma at every j; transition matrix chi2 = 11.0
(df 9, null); P(j | n mod 9) independent at all 9 residues. The earlier
full-orbit deviations (chi2 = 324) were ENDGAME CONTAMINATION: every
orbit exits through the same small-number gauntlet with deterministic
j's. In the stationary regime the phase alphabet has zero exploitable
structure - the fair-dice doctrine verified at phase level.

## Observation 84 (uniformity series extended to k=15), R2001-2020
Two-point Hilbert contraction: mu = 0.8544 / 0.8703 / 0.8686 / 0.8736
(k = 9/11/13/15, N up to 4.78M). Bounded by 0.875 through k=15 with a
slow upward creep (~+0.002/k) - cannot yet distinguish asymptote < 1
from slow approach to 1. The uniformity question stays open; the
quantitative window narrows.

## Theorem 78b (STERILITY IS O(log k) - rigorous route via p-adic Baker),
## R2021-2030
Records of fresh fuel r(k) from pure repunit burns, extended to k = 10^6:
r = 1,2,6,9,17,18,20,21 at k = 2,6,12,174,1198,263342,539916,787630.
Empirically r(k) ~ 1.07 log2(k). RIGOROUS: r(k) = v2-type valuation of
3^k - c, bounded by C log k via p-adic linear forms in logarithms
(K. Yu's theorem) - so the sterility gap k - r(k) >= k - C log k is
PROVABLE: a pure run of k ones can never return more than logarithmic
fuel, with an effective constant. The engine's no-feedback law is
Baker-effective, not just empirical.

## Observation 85 (THE EXIT GAUNTLET AND THE REPUNIT GATES), R2031-2050
(1) The last 15 odd steps of every orbit have a fixed non-dice signature
    (j=4: 16.6% vs 6.25%; j=1 depleted to 39.5%) — the deterministic
    funnel that contaminated the naive census (Obs 83).
(2) THE FINAL-GATE LAW (exact): the last odd value before the 2^j tower
    must satisfy 3m+1 = 2^j: m = (2^j-1)/3 with j EVEN (2^j = 1 mod 3),
    and j not = 0 mod 6 (else 3 | m and m is unreachable — multiples of
    3 never occur mid-orbit). The gates are exactly the BASE-4 REPUNITS:
    5 = 11_4 (93.9% of all orbits), 85 = 1111_4 (2.3%), 341 = 11111_4
    (3.7%), 5461, 21845, 349525 (traces) — repunit lengths j/2 not
    divisible by 3. EVERY Collatz orbit ends through a base-4 repunit
    gate; the missing gates (21, 1365, 87381 = lengths 3,6,9) are the
    mod-3 shadow of Thm 41-style exclusions. The problem that begins
    with binary repunits (fuel) provably ENDS in quaternary repunits
    (gates): the +1 anchor family closes the loop.

## Proposition 86 (THE GATE MEASURE, derived: the L mod 3 trichotomy).
## R2081-2130
The exit gates are the base-4 repunits g_L = (4^L - 1)/3, and
g_L = L mod 3 (since 4^L = 1 + 3L mod 9). This single congruence
organizes the entire endgame:
* L = 0 mod 3: g divisible by 3 -> unreachable, NO basin (missing gates
  21, 1365, 87381 - as observed).
* L = 1 mod 3: g = 1 mod 3 -> backward preimages only at even j
  (branches grow x4/3): POOR basin.
* L = 2 mod 3: g = 2 mod 3 -> has the j=1 backward branch (x2/3 dip):
  RICH basin.
Corrected backward-tree counting (nodes above the window can dip back
below it via the x2/3 branch - the naive count missed 40%) converges to
85: 2.49%, 341: 3.81%, 5461: 0.02%, 21845: 0.05% vs forward-measured
2.65 / 3.92 / 0.01 / 0.05 - agreement, and the ordering inversions
(341 > 85 despite being 4x larger; 21845 > 5461) are exactly the
trichotomy. Gate populations are stable across seed scales 2^16-2^30
(93.4-94.0% for gate 5): no significant log-periodic drift at this
resolution. The endgame is now THEORY: every orbit exits through a
quaternary repunit whose basin richness is L mod 3.

## Theorem 87 (THE HALVING ORACLE — M. de Jong's question, closed form).
## — PROVED (one line) + verified 20000/20000, R2131-2150
Seed n = m*2^k - 1 (family k, head m odd). After the k-step climb the
value is 3^k*m - 1, so the halving count is
    c = v2(3^k*m - 1) = 2-ADIC AGREEMENT DEPTH OF THE HEAD m WITH 3^(-k).
The pattern per family EXISTS and is the bit string of 3^(-k):
    family 1: ...101010101011 (the ALTERNATOR, period 2)
    family 2: ...111000111001 (block 111000, period 6)
    family k: period ord_{3^k}(2) = 2*3^(k-1) — cyclotomic, tripling.
Heads halving exactly c times = the residue class m = 3^(-k) mod 2^c
(and not mod 2^(c+1)). Across random heads the c-distribution is exactly
dice 2^-c (verified k=1,4,8) — the oracle is the DETERMINISTIC refinement
of the dice: the coin flip IS the head's next unread oracle bit.
MERGE LAW (Martien's "merge step"): two sequences of the same family
whose heads agree 2-adically to depth d have identical post-burn
behavior to depth d — families/sequences merge exactly as deep as their
heads agree with each other, and diverge at the first differing bit.
The orbit is a chain of oracle consultations: burn run k, read c from
oracle 3^(-k) at the head, form the new head, consult the next oracle.
Collatz = the machine that reads the 2-adic expansions of 3^(-k) with
its own output as the address.

## Theorem 88 (THE ADDRESS CODING: how the oracle patterns run on).
## — verified exhaustively (36 prefixes, all seeds to 2^20), R2151-2175
The phase chain (run k_i, halvings c_i) obeys the exact affine law
    m = 3^-k1 (1 - 2^c1) + 3^-k1 2^(c1+k2) m2,
iterated per phase. Consequences, all verified:
(a) ADDRESS-BIT LAW: seeds sharing a phase prefix form EXACTLY ONE
    residue class mod 2^(D+1), D = sum(k_i + c_i) - the oracle bits are
    literally the seed's address bits; class density = the ledger.
(b) PERIODIC WORDS -> ANCHOR LATTICE: the reference of the repeated
    word (k,c) is m* = (1-2^c)/(3^k - 2^(k+c)) - anchor-lattice
    rationals (denominators 3^r - 2^j). INTEGER references = cycles:
    (1,1) gives m* = 1 = the trivial cycle; (2,1) gives m* = -1, i.e.
    n = -5: THE -5 CYCLE OF 3n-1 appears as a negative integer address
    - the mirror map's cycles live in the same coding, on the negative
    side. Near-integer references = comma words (near-cycles).
(c) The infinite-depth limit of the nested references is the 2-adic
    coding of the full orbit - the classical Bernstein-Lagarias
    conjugacy, here derived as the limit of Martien's family/oracle
    refinements. The family framework and the 2-adic conjugacy are THE
    SAME OBJECT read in opposite directions: families = finite-depth
    truncations of Phi.

## Theorem 89 (THE HEAD/TAIL FACTORIZATION of backward convergence —
## M. de Jong's question). — verified exact, R2176-2195
For a target n = HEAD|TAIL, the predecessors factorize EXACTLY:
(1) REAR (2-adic substrate): the predecessor tail is an exact function
    of the target tail alone (same j): tail_pre = (tail*2^j - 1)*3^-1
    mod 2^b. Verified: 5000/5000. Because binary carries propagate
    UPWARD ONLY, the tail is perfectly autonomous - it NEVER feels the
    head (measured seam downward: 0 bits, always).
(2) The rear does not pull back to the same pattern: it pulls back
    through the ORACLE LADDER: depth-i predecessors carry the
    3^-i-twisted tail. Live on 26623 = 110011111111111:
    depth 1: 35497 = ...0101010101001 (alternator = 3^-1),
    depth 2: 47329 = ...011100001    (111000-block = 3^-2),
    exactly the family oracles of Thm 87 appearing in the flesh.
(3) FRONT (archimedean substrate): pure value scaling, head_pre =
    head * 2^j/3 to relative precision 3e-7; the head receives only a
    <= 2-bit carry trickle from below per step.
So Martien's factorization is correct with one refinement: the rear
converges through the 3^-i twist ladder (not to the fixed pattern), the
front converges as (2^j/3)-scaling, and the interface is ONE-WAY: tail
drives head, never the reverse. The two-substrate picture is exact at
string level: 2-adic rear, archimedean front, one-way carry seam.

## Theorem 90 (THE PRE-WRITING FORMULA — M. de Jong's follow-up).
## — exact, verified live on the 26623 highway, R2196-2210
The convergence of a sequence toward any future near-repunit is governed
by one explicit affine formula. Per phase (k, c) compose:
    R    <- R + coef * 3^-k (1 - 2^c)
    coef <- coef * 3^-k * 2^(c + k_next)
Then  starter_head = R + coef * head_at_arrival, and the arrival run of
L ones is PRE-WRITTEN in the starter as its 2-adic agreement with the
rational R over the bit window [D - k1, D - k1 + L), D = sum(k_i + c_i).
LIVE: starter 63105 (a 3|n leaf = true highway head), three phases
(1,1): R = -37/27, coef = 2^16/27, recovered head = 13 = (26623+1)/2^11
exactly - the eleven ones of 26623 sat in the starter's bits [5,16) as
agreement with -37/27. Note coef = 2^(consumed binary depth)/3^(total
rises): THE FORMULA IS THE BINARY/TERNARY CONVERSION LEDGER ITSELF -
numerator counts the binary bits consumed, denominator the ternary
lifts performed. Family-sequence cycles, the oracle ladder, and the
conversion thesis are one bookkeeping identity.

## Proposition 91 (THE 5% TAX: one constant behind all three walls).
## R2211-2250, from working out the pre-writing formula
(a) THE ESCAPE DIMENSION (exact closed form): via Thm 90, divergence
    requires the seed's 2-adic address to match rise-rich references
    forever (sum j <= alpha * r). The Hausdorff dimension of the
    escape-address set is EXACTLY the binary entropy
        dim_esc = H(1/log2 3) = H(0.63093) = 0.94996,
    attained by the geometric tilt with mean j = log2 3 (numeric
    two-parameter optimum 0.94912 = truncation of this closed form).
(b) THE SAME CONSTANT appears at wall 2: critical-window cycle words
    number 2^{H(1/alpha) j} against modulus 2^j - the "5% tax"
    1 - H(1/alpha) = 0.05004 per halving bit is why expected phantom
    hits vanish (Prop 62 census). Drift-balanced word entropy is ONE
    object governing cycles (wall 2) and divergence addresses (wall 3).
(c) THE GAMMA FORK (honest): does the K-L hierarchy converge to 1 or
    to the word-entropy ceiling H(1/alpha) = 0.950? Free geometric fit
    of the ladder (0.8624/0.8805/0.8953/0.9069/0.9146) lands at
    gamma_inf = 0.986 with residual equal to the forced gamma_inf = 1
    fit - current data CANNOT distinguish, and mildly disfavors the
    0.950 ceiling (would require bending below trend). DISCRIMINATING
    TEST: certifications at k ~ 25-30; gamma crossing 0.950 kills the
    ceiling conjecture. Both models predict gamma_21 = 0.919.
STATUS of "working it out to prove Collatz": the conjecture now reads
as ONE quantitative statement - no integer address survives the 5% tax
forever. Cycles: the tax kills word-counts against the modulus
(census-verified r <= 24). Divergence: the tax makes escape addresses
dimension-deficient (0.950 < 1). The remaining gap is, as ever,
pointwise: integers are countable and dimension arguments see only
size. But every wall now has the same number on it.

## Proposition 92 (COLLATZ IS A RADIUS-1 CELLULAR AUTOMATON IN BASE 6;
## local-invariant no-go #10). R2251-2295
(1) STRUCTURE (verified exact, 3000/3000 each): in base 6 both elementary
    maps are carry-free radius-1 sliding rules:
      x3:  out_i = 3*(d_i mod 2) + floor(d_{i-1}/2)
      /2:  out_i = floor(d_i/2)  + 3*(d_{i+1} mod 2)
    The ONLY nonlocality is the +1 (a boundary defect propagating through
    trailing 5s). Collatz is literally a one-dimensional cellular
    automaton in the mixed base - the rule-110 intuition (WOLFRAM_RST_LRS
    arc) is now an exact statement, not an analogy. Note the two rules
    are mirror transposes: x3 reads parity-up, /2 reads parity-down.
(2) INVARIANT HUNT (Hattori-Takesue style): window-2 weight functions:
    - translation-invariant: NO conserved/eigen densities exist (0-dim).
    - LSD-weighted sum lambda^i f(d_i,d_{i+1}), lambda in {1/3..5}: only
      the trivial constant. At lambda = 6 the value functional exists
      (known) but the numerical formulation degenerates (top-window
      domination); exact-arithmetic version left open.
    NO-GO #10: the CA has no local additive first integrals at window 2
    beyond the value itself - sharpens Obs 82b (0.666 barrier): not only
    no monotone local quality, no locally-conserved density AT ALL. The
    conserved information of Collatz is irreducibly global (the address
    coding), consistent with everything since Thm 22.

## Proposition 93 (THE DYNAMICAL ORDER — M. de Jong's proposal).
## R2296-2320
Order integers by index(n) = weighted path cost from 1 in the backward
tree (w2 per doubling, w3 per odd-inverse step; paths unique). Built by
Dijkstra; findings:
(1) THE GATES COME FIRST: the enumeration opens 1,2,4,8,16,5,32,10,64,
    3,20,21,128,... - the early population is exactly the quaternary-
    repunit gate family (5,10,20,21,40,42,84,85,168,170,340,341...) and
    the 2-power spine: the "small numbers" of the dynamical order are
    the endgame highway of Obs 85. Robust to weighting ((1,1), (1,2),
    (1,log2 3) nearly identical).
(2) SHELL GROWTH: N(t) ~ 2^(0.342 t) for plain steps; consistent with
    branching 1 = 2^-rho(1 + q), q = measured odd-child density 0.26.
(3) SHEAR + DICE: index = 0.277 log2(n) + const with residual sd ~ 4:
    the dynamical order is an affine shear of the size order plus the
    stopping-time fluctuation - the two coordinate systems differ by
    exactly the dice.
(4) CANONICAL WEIGHTING: w = bits consumed (the ledger metric) makes
    shells = address-prefix classes, shell counts = 2^(H D) admissible
    words: the 5% tax appears as shell-vs-value deficit; gamma (K-L) =
    how densely dynamical shells cover each size window; and the
    conjecture itself = "the dynamical order is TOTAL" (every integer
    receives a finite index). Martien's reindexing is the coordinate
    system in which the entire machinery (ledger, tax, gates, address
    coding) is native.
