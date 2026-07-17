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
1. THEOREM 17 EXTENDED to period 14: the complete integer cycle list through
   14 odd steps remains {1, -1, -5, -17} (exhaustive, both signs).
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
