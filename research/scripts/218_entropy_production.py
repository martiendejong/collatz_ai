"""
218_entropy_production.py
=========================
Methode 9: Entropieproductie als bewijs-instrument.

De Collatz-map T is niet tijdsomkeerbaar: elk oneven getal n heeft
meerdere voorgangers (de inverse T^{-1}(n) bevat meerdere getallen).

Shannon-entropie van de inverse map:
  H_inv(n) = -Σ_m P(T(m)=n) log P(T(m)=n)

Als H_inv(n) > 0 voor alle n: de map produceert informatie per stap
-> thermodynamisch argument tegen cycli.

Binnen de K-L-structuur: voor residuen mod 3^k:
  - Elke knoop i heeft een verzameling van voorgangers in de boom
  - De entropieproductie = log(|{j : T(j) = i}|) gemiddeld

We meten:
  (A) De verdelingsbreedte van het aantal voorgangers per knoop
  (B) De gemiddelde Shannon-entropie van de inverse map
  (C) Of de entropieproductie uniform begrensd is weg van nul

Verband met cyclus-vrijheid: als H_inv(n) >= delta > 0 voor alle n,
dan is de "thermodynamische weerstand" tegen terugkeer naar n
minstens delta bits per stap. Voor een cyclus van lengte k geldt dan:
  totale entropieproductie >= k * delta > 0
maar een cyclus heeft nul netto entropieproductie (per definitie
terugkeer) -> tegenspraak voor delta > 0 uniform.

Dit is geen volledig bewijs, maar geeft de thermodynamische intuïtie
een kwantitatieve basis.
"""
import numpy as np
from math import log2, log
from collections import defaultdict

ALPHA = log2(3.0)
LAM   = 1.70
A     = LAM ** -2.0
B1    = LAM ** (ALPHA - 2.0)
B3    = LAM ** (ALPHA - 1.0)


def perron(k, n_iter=300):
    N  = 3 ** (k - 1)
    i  = np.arange(N, dtype=np.int64)
    T4 = (4 * i + 2) % N
    s, r = np.divmod(i, 3)
    Nl = N // 3
    m0, m2 = (r == 0), (r == 2)
    R1 = (4 * s) % Nl
    R3 = (2 * s + 1) % Nl
    v = np.ones(N, dtype=np.float64)
    for _ in range(n_iter):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        w  = A * v[T4]
        w[m2] += B3 * cb[R3[m2]]
        w[m0] += B1 * cb[R1[m0]]
        v = w / w.max()
    v /= v.sum()
    return v, T4, R1, R3, m0, m2, Nl, N


def count_predecessors(k):
    """
    Voor elke knoop i: tel het aantal voorgangers in de K-L-boom.
    Een voorganger j -> i bestaat als T4(j) = i of R3(j) = i of R1(j) = i.
    """
    v, T4, R1, R3, m0, m2, Nl, N = perron(k)
    i_all = np.arange(N, dtype=np.int64)

    # Voorwaartse verbindingen:
    # Walk: j -> T4(j)  voor alle j
    # Lift type 2: j -> R3(j) + sel[R3(j)]*Nl  (maar dit is de gelinieariseerde versie)
    # In de eigenlijke boom: de T4-map is de enige structurele verbinding
    # die op het niveau van de volledige residu-structuur werkt.

    # Eenvoudiger: tel hoeveel j's naar elke i mappen via T4
    pred_count = np.zeros(N, dtype=np.int64)
    for j in i_all:
        pred_count[int(T4[j])] += 1

    # Aanvullende voorgangers via de lift-operatoren
    # (voor elke type-2 knoop j: ook R3(j) heeft j als voorganger)
    lift2_targets = (R3 + 0 * Nl)[m2]  # vereenvoudigd: eerste sel=0
    for tgt in lift2_targets:
        if 0 <= tgt < N:
            pred_count[tgt] += 1

    lift0_targets = R1[m0]
    for tgt in lift0_targets:
        if 0 <= tgt < N:
            pred_count[tgt] += 1

    return pred_count, v


def entropy_production_per_node(pred_count, mu):
    """
    H_inv(i) = log(pred_count[i])  [als alle voorgangers gelijkwaardig]
    Gewogen gemiddelde: E_mu[H_inv] = Σ_i mu_i * log(pred_count[i])
    """
    eps = 1e-10
    h = np.where(pred_count > 0, np.log(pred_count.astype(float) + eps), 0.0)
    return h, float(np.dot(h, mu))


# Collatz reverse map op echte getallen (voor vergelijking)
def collatz_predecessors(n, max_pred=20):
    """Vind voorgangers van n in de echte Collatz-stroom."""
    preds = []
    # Type 1: n = m/2 -> m = 2n (als m niet door 3 deelbaar is via T(m)=n)
    m = 2 * n
    preds.append(m)
    m = 4 * n
    preds.append(m)
    # Type 2: n = (3m+1)/2^k -> 3m = 2^k*n - 1 -> m = (2^k*n-1)/3
    for k_pow in range(1, 10):
        val = (2**k_pow * n - 1)
        if val > 0 and val % 3 == 0:
            m = val // 3
            if m % 2 == 1:  # m oneven
                preds.append(m)
    return preds[:max_pred]


print(f"Methode 9: Entropieproductie  (lam={LAM})", flush=True)
print("=" * 65, flush=True)

print("\n(A) K-L booms-structuur: voorgangersverdeling", flush=True)
for k in (11, 12, 13):
    pred_count, mu = count_predecessors(k)
    h, E_h = entropy_production_per_node(pred_count, mu)
    print(f"  k={k}  N={3**(k-1):,}", flush=True)
    print(f"    pred_count: min={pred_count.min()}  max={pred_count.max()}  "
          f"mean={pred_count.mean():.3f}", flush=True)
    print(f"    E_mu[log(preds)] = {E_h:.5f}  (entropieproductie/stap)", flush=True)
    print(f"    P(pred_count=1) = {(pred_count==1).mean():.4f}  "
          f"P(>=3) = {(pred_count>=3).mean():.4f}", flush=True)
    print(f"    H_inv: min={h.min():.4f}  max={h.max():.4f}  "
          f"std={h.std():.4f}", flush=True)

print("\n(B) Echte Collatz: voorgangersverdeling voor grote n", flush=True)
# Steekproef van oneven getallen
rng  = np.random.default_rng(42)
ns   = rng.integers(1, 10**6, size=2000) * 2 + 1  # oneven
n_preds = [len(collatz_predecessors(int(n))) for n in ns]
h_vals  = [np.log(max(p, 1)) for p in n_preds]
print(f"  E[#voorgangers] = {np.mean(n_preds):.3f}", flush=True)
print(f"  E[log(#voorgangers)] = {np.mean(h_vals):.5f}", flush=True)
print(f"  Min voorgangers = {min(n_preds)}  Max = {max(n_preds)}", flush=True)

print("\n(C) Thermodynamisch argument:", flush=True)
k_mid = 12
pred_count, mu = count_predecessors(k_mid)
min_preds = pred_count.min()
E_h = float(np.dot(np.log(pred_count.astype(float) + 1e-10), mu))
print(f"  Minimaal voorgangers per knoop (k={k_mid}): {min_preds}", flush=True)
print(f"  Gemiddelde entropieproductie: {E_h:.5f} nat/stap", flush=True)
print(f"  = {E_h/log(2):.5f} bits/stap", flush=True)
print(f"  Implicatie: elke cyclus van lengte L produceert >= "
      f"{min_preds} * L bits entropie", flush=True)
print(f"  Maar cyclus heeft NETTO 0 entropieproductie -> contradictie als min>=2?",
      flush=True)
print(f"  (Nee: het GEMIDDELDE is positief, niet elke stap; formele kloof blijft.)",
      flush=True)

print("\ndone", flush=True)
