"""
220_perturbation_theory.py
==========================
Eerste-orde storingstheorie voor de K-L-operator.

Bij eps=0: L_0 = A * T4 (zuivere walk-op-4, geen liften)
Eigenvector v^0 = (1,1,...,1) (uniform), eigenwaarde rho_0 = A = lam^{-2}

Bij eps>0: L = L_0 + eps * L_1
  L_1 = liftoperator: [L_1 v](i) = B3*v[R3(i)] als r=2,
                                     B1*v[R1(i)] als r=0, 0 anders
  (exact: L_1 werkt op het minimum cb, maar bij eps->0 is cb=min(v^0)=1 overal)

Eerste-orde storingstheorie (niet-degenereerte eigenwaarde):
  rho(eps) = rho_0 + eps * <u^0, L_1 v^0> / <u^0, v^0>  + O(eps^2)
  v(eps)   = v^0   + eps * V_1 + O(eps^2)

waarbij V_1 de oplossing is van:
  (L_0 - rho_0 * I) V_1 = (rho_1 * I - L_1) v^0
  met rho_1 = <u^0, L_1 v^0> / <u^0, v^0>

ANALYTISCH DOEL:
  Toon V_1 != 0 aan -> d(CV)/d(eps)|_{eps=0} = std(V_1) > 0 -> via lineariteit:
  CV(eps=1) >= d(CV)/d(eps)|_{eps=0} * 1 > 0

BEREKENING:
  L_0 v = A * v[T4]  (T4 is een PERMUTATIE van Z/N*Z, want 4 is inverteerbaar mod 3)
  Dus L_0 heeft ALLE EIGENWAARDEN gelijk aan A (want T4 is een bijekties permutatie
  van de coordinaten en A is een scalaire factor).

  Wacht: dat klopt niet. T4: i -> (4i+2) mod N is een PERMUTATIE, dus L_0 v = A * P v
  waarbij P de T4-permutatiematrix is.
  De eigenwaarden van P zijn de N-de en lagere eenheidswortels.
  L_0 heeft eigenwaarden A * zeta waarbij zeta over de eigenwaarden van P loopt.
  De Perron-eigenwaarde van L_0 is A (voor de uniforme eigenvector, en T4 heeft
  eigenwaarde 1 voor de uniforme vector omdat T4 een bijectie is).

INVERSEBEREKENING:
  (L_0 - A*I) V_1 = -(L_1 - rho_1*I) v^0 = -L_1 v^0 + rho_1 * ones

  L_0 - A*I = A*(P - I)

  Vereist: (P - I) V_1 = -(1/A) * (L_1 v^0 - rho_1 * ones)

  Omdat P een permutatie is, heeft P - I een nulruimte = span(ones)
  (Pones = ones omdat P een bijectie is). Dus (P-I)^{-1} bestaat op de ruimte
  loodrecht op ones.

  De RHS -(1/A)*(L_1 ones - rho_1 * ones) heeft gemiddelde nul als we rho_1
  kiezen gelijk aan het gemiddelde van L_1 ones.

  V_1 = UNIQUE oplossing in de orthogonale complement van ones.
"""
import numpy as np
from math import log2

ALPHA = log2(3.0)
LAM   = 1.70
A_coef = LAM ** -2.0
B1_coef = LAM ** (ALPHA - 2.0)
B3_coef = LAM ** (ALPHA - 1.0)


def build_T4(k):
    N = 3 ** (k - 1)
    i = np.arange(N, dtype=np.int64)
    T4 = (4 * i + 2) % N
    return T4, N


def lift_action(v, k):
    """
    L_1 v bij eps=0: cb = min(v0,v1,v2) = v (uniform bij eps=0),
    maar voor willekeurige v: cb(s) = min(v[s], v[s+Nl], v[s+2*Nl])
    """
    N = v.size
    i = np.arange(N, dtype=np.int64)
    s, r = np.divmod(i, 3)
    Nl = N // 3
    m0, m2 = (r == 0), (r == 2)
    R1 = (4 * s) % Nl
    R3 = (2 * s + 1) % Nl
    # Bij eps=0 is v uniform, dus cb = v[:Nl] = v[Nl:] = v[2*Nl:] = 1
    # We berekenen L_1 v = B3*cb[R3] als r=2, B1*cb[R1] als r=0
    cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
    w  = np.zeros_like(v)
    w[m2] = B3_coef * cb[R3[m2]]
    w[m0] = B1_coef * cb[R1[m0]]
    return w


def solve_perturbation(k):
    """
    Berekening van de eerste-orde correctie V_1.

    (L_0 - A*I) V_1 = -(L_1 - rho_1*I) v^0

    L_0 = A * P_T4  (permutatiematrix voor T4)
    Oplossing via Fourier op de permutatie-groep:
      V_1[i] = som over cycli van T4 van de RHS-waarden
    """
    T4, N = build_T4(k)
    v0 = np.ones(N, dtype=np.float64)

    # L_1 v^0 (v^0 = uniform = ones)
    L1_v0 = lift_action(v0, k)
    rho_1 = float(L1_v0.mean())
    print(f"  k={k}:  rho_1 = {rho_1:.7f}  (A={A_coef:.7f}  rho_echt~{A_coef+rho_1:.7f})")

    # RHS = -(L_1 v^0 - rho_1 * ones) = rho_1 - L_1 v^0
    rhs = rho_1 - L1_v0  # heeft gemiddelde 0 per constructie

    # Oplossing van A*(P - I) V_1 = rhs via  (P - I) V_1 = rhs / A
    # P - I: itereer de permutatie T4 om de cycli te vinden,
    # dan V_1[i] = -(rhs[i] + rhs[T4[i]] + ... + rhs[T4^{l-1}[i]]) / (A * l)
    # voor elke cyclus van lengte l.

    rhs_scaled = rhs / A_coef  # = (rho_1 - L_1 v^0) / A

    # Bereken cycli van T4
    visited = np.zeros(N, dtype=bool)
    V1 = np.zeros(N, dtype=np.float64)

    cycle_lengths = []
    for start in range(N):
        if visited[start]:
            continue
        # Volg de cyclus
        cyc = []
        cur = start
        while not visited[cur]:
            cyc.append(cur)
            visited[cur] = True
            cur = int(T4[cur])
        L_cyc = len(cyc)
        cycle_lengths.append(L_cyc)
        # Cyclische som: (P - I) V_1 = rhs/A  binnen de cyclus
        # Als de cyclus [i0, i1, ..., i_{L-1}] is (i_{j+1} = T4(i_j)):
        # V_1[i_{j+1}] - V_1[i_j] = rhs_scaled[i_j]
        # Oplossing: V_1[i_j] = V_1[i_0] + sum_{l=0}^{j-1} rhs_scaled[i_l]
        # Consistentie: sum_{l=0}^{L-1} rhs_scaled[i_l] = 0 (vereist!)
        cyc_arr = np.array(cyc)
        cyc_rhs = rhs_scaled[cyc_arr]
        cyc_sum = float(cyc_rhs.sum())
        if abs(cyc_sum) > 1e-8:
            # Inconsistente cyclus -> orthogonaalprojectie correctie
            cyc_rhs -= cyc_sum / L_cyc
        # Cumulatieve som geeft V_1
        V1_cyc = np.cumsum(cyc_rhs)
        V1_cyc -= V1_cyc.mean()
        V1[cyc_arr] = V1_cyc

    V1 -= V1.mean()
    cv_V1 = float(np.std(V1))  # d(CV)/d(eps) bij eps=0
    print(f"  std(V_1) = {cv_V1:.6f}  (= d(CV)/d(eps) bij eps=0)")
    print(f"  max|V_1| = {np.abs(V1).max():.6f}")

    # Cycli-statistiek
    cycle_lengths = np.array(cycle_lengths)
    print(f"  T4-cycli: {len(cycle_lengths)} cycli, "
          f"lengte min={cycle_lengths.min()} max={cycle_lengths.max()} "
          f"mean={cycle_lengths.mean():.2f}")
    # Langste cyclus bepaalt de structuur
    uniq, cnt = np.unique(cycle_lengths, return_counts=True)
    print(f"  Cycli-verdeling: {dict(zip(uniq[:5].tolist(), cnt[:5].tolist()))}")

    # Verificatie: corr(V1, v_eps-v0) bij kleine eps
    return V1, cv_V1, rho_1


def verify_with_perron(k, eps=0.05):
    """Vergelijk V_1 met (v(eps) - v^0) / eps."""
    N  = 3 ** (k - 1)
    i  = np.arange(N, dtype=np.int64)
    T4 = (4 * i + 2) % N
    s, r = np.divmod(i, 3)
    Nl = N // 3
    m0, m2 = (r == 0), (r == 2)
    R1 = (4 * s) % Nl
    R3 = (2 * s + 1) % Nl
    B1 = eps * B1_coef
    B3 = eps * B3_coef

    v = np.ones(N, dtype=np.float64)
    for _ in range(300):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        w  = A_coef * v[T4]
        w[m2] += B3 * cb[R3[m2]]
        w[m0] += B1 * cb[R1[m0]]
        v = w / w.max()
    v -= v.mean(); v /= eps
    return v


print(f"220: Eerste-orde storingstheorie  (lam={LAM})")
print("=" * 65)
print()

for k in (8, 9, 10, 11, 12):
    V1, cv_V1, rho_1 = solve_perturbation(k)

    # Verificatie bij k<=10
    if k <= 10:
        eps_test = 0.02
        v_num = verify_with_perron(k, eps_test)
        corr_V1 = float(np.corrcoef(V1, v_num)[0, 1])
        print(f"  corr(V_1, (v(eps)-v0)/eps) = {corr_V1:.5f}  (moet ~1 zijn)")
    print()

print("Samenvatting d(CV)/d(eps) bij eps=0:")
print("  std(V_1) is de linearisatie-coefficient.")
print("  Als std(V_1) > 0 voor alle k, en CV(eps) lineair in eps,")
print("  dan CV(eps=1) = std(V_1) * 1 > 0 -> G.")
print()
print("done")
