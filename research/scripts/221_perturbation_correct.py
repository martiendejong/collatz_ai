"""
221_perturbation_correct.py
===========================
Gecorrigeerde eerste-orde storingstheorie voor de K-L-operator.

Het probleem in script 220: bij eps=0 is v=1 (uniform), maar
min(1,1,1) is NIET differentieerbaar — alle drie componenten zijn gelijk.
De correcte linearisatie bij de SYMMETRIEPUNT is:
  d/d(delta_v) min(1+eps*a, 1+eps*b, 1+eps*c)|_{eps=0} = min(a,b,c)
Maar als a=b=c (uniform perturbatie), geeft de min de laagste.

Oplossing: gebruik de SUBGRADIËNT op het symmetriepunt.
De convexe subgradiënt van min(x,y,z) bij x=y=z=1 is
de convexe hul van alle gradients: de gewichten (w1,w2,w3) met
w1+w2+w3=1, alle wi>=0, mogen willekeurig verdeeld zijn.

In de context van storingstheorie: als we v = 1 + eps*u perturb,
dan is cb(s) = min(v[s], v[s+Nl], v[s+2Nl])
            = min(1+eps*u[s], 1+eps*u[s+Nl], 1+eps*u[s+2Nl])
            = 1 + eps * min(u[s], u[s+Nl], u[s+2Nl])  + O(eps^2)

Dus de gecorrigeerde linearisatie van de liftterm is:
  L_1^correct v = B3 * min(v[s], v[s+Nl], v[s+2Nl]) (voor r=2)
                + B1 * min(v[s], v[s+Nl], v[s+2Nl]) (voor r=0)

Dit is NIET lineair in v! Maar we kunnen het benaderen via:

METHODE A: Directe numerieke berekening van d(CV)/d(eps) via
  V_1 = lim_{eps->0} (v(eps) - v(0)) / eps

  waarbij v(eps) de Perron-eigenvector is bij liftsterkte eps.
  Dit vermijdt het linearisatieprobleem volledig.

METHODE B: Gebruik dat bij eps=0 v=uniform, dus min(u[s],u[s+Nl],u[s+2Nl])
  exact berekend kan worden voor elke perturbatierichting u.
  De operator L_1^lin: u -> B3*min(u[s],u[s+Nl],u[s+2Nl]) is
  een piecewise lineaire operator op perturbaties.
"""
import numpy as np
from math import log2

ALPHA = log2(3.0)
LAM   = 1.70
A_coef  = LAM ** -2.0
B1_coef = LAM ** (ALPHA - 2.0)
B3_coef = LAM ** (ALPHA - 1.0)


def perron_eps(k, eps, n_iter=300):
    """Perron-vector bij liftsterkte eps."""
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
    for _ in range(n_iter):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        w  = A_coef * v[T4]
        if eps > 0:
            w[m2] += B3 * cb[R3[m2]]
            w[m0] += B1 * cb[R1[m0]]
        v = w / w.max()
    v /= v.mean()
    return v


def numerical_derivative(k, eps_vals=(0.005, 0.010, 0.020, 0.050)):
    """Numerieke afgeleiden van CV(eps) bij eps->0."""
    v0 = perron_eps(k, 0.0)
    print(f"  k={k}  N={3**(k-1)}")
    for eps in eps_vals:
        v_eps = perron_eps(k, eps)
        # Directe perturbatievector
        dv = (v_eps - v0) / eps
        dv -= dv.mean()
        std_dv = float(np.std(dv))
        corr = float(np.corrcoef(v_eps, v0)[0, 1])
        # d(CV)/d(eps) bij dit eps
        dcv_deps = float(np.std(v_eps) / eps)  # CV(eps)/eps (voor CV(0)=0)
        print(f"    eps={eps:.3f}: std(dv/eps)={std_dv:.4f}  "
              f"CV/eps={dcv_deps:.4f}  corr(v_eps, v0)={corr:.4f}")
    print()


def correct_V1(k, n_iter_power=300):
    """
    Bereken de echte eerste-orde correctie via de gecorrigeerde lineaire operator.

    (L_0 - A*I) V_1 = -(L_1^lin - rho_1*I) ones

    waarbij L_1^lin u = B3 * min(u[s], u[s+Nl], u[s+2Nl]) voor r=2
                      + B1 * min(u[s], u[s+Nl], u[s+2Nl]) voor r=0

    Maar dit is NON-LINEAIR in u! De correct gelinieariseerde versie
    rond v^0 = ones is:
      L_1^lin u[i] = B3 * u[R3[i] + argmin_j(ones) * Nl] voor r=2
      = B3 * u[R3[i]] als alle ties breken naar de EERSTE component

    Maar bij uniform v^0=ones is argmin ALLE drie gelijk.
    We kiezen de LAGRANGIAN linearisatie: voor een uniforme perturbatie u,
    de min-linearisatie geeft GEMIDDELD (u[s]+u[s+Nl]+u[s+2Nl])/3.

    Dit zijn twee uitersten:
    (I) min-linearisatie: L_1^lin u = B3 * min(u[s], u[s+Nl], u[s+2Nl])
    (II) gem-linearisatie: L_1^lin u = B3 * (u[s]+u[s+Nl]+u[s+2Nl])/3
    (III) directe numerieke V_1 via finite difference
    """
    N  = 3 ** (k - 1)
    i  = np.arange(N, dtype=np.int64)
    T4 = (4 * i + 2) % N
    s, r = np.divmod(i, 3)
    Nl = N // 3
    m0, m2 = (r == 0), (r == 2)
    R1 = (4 * s) % Nl
    R3 = (2 * s + 1) % Nl

    # Methode III: finite difference V_1 = (v(eps) - 1) / eps voor kleine eps
    eps_small = 0.005
    v_eps = perron_eps(k, eps_small, n_iter=500)
    V1_fd = (v_eps - 1.0) / eps_small
    V1_fd -= V1_fd.mean()
    std_fd = float(np.std(V1_fd))

    # Methode I: min-linearisatie
    def L1_min(u):
        """L_1^lin u met min-linearisatie."""
        cb_u = np.minimum(np.minimum(u[:Nl], u[Nl:2*Nl]), u[2*Nl:])
        w = np.zeros_like(u)
        w[m2] = B3_coef * cb_u[R3[m2]]
        w[m0] = B1_coef * cb_u[R1[m0]]
        return w

    # Methode II: gemiddelde-linearisatie
    def L1_avg(u):
        """L_1^lin u met gemiddelde-linearisatie."""
        avg_u = (u[:Nl] + u[Nl:2*Nl] + u[2*Nl:]) / 3.0
        w = np.zeros_like(u)
        w[m2] = B3_coef * avg_u[R3[m2]]
        w[m0] = B1_coef * avg_u[R1[m0]]
        return w

    # rho_1 (beide methoden geven hetzelfde voor v^0 = ones!)
    ones = np.ones(N, dtype=np.float64)
    rho_1 = float(L1_min(ones).mean())  # = float(L1_avg(ones).mean()) = (B3+B1)/3

    # Oplossing (L_0 - A*I) V_1 = -(L_1 - rho_1 I) ones = rho_1 * ones - L_1 ones
    # Via T4-cyclus methode (uit script 220):
    def solve_via_cycle(rhs_vec):
        """Los (A*(P_T4 - I)) V op voor RHS = rhs_vec via cyclisch cumulatief som."""
        rhs_scaled = rhs_vec / A_coef
        visited = np.zeros(N, dtype=bool)
        V = np.zeros(N, dtype=np.float64)
        for start in range(N):
            if visited[start]: continue
            cyc = []
            cur = start
            while not visited[cur]:
                cyc.append(cur); visited[cur] = True; cur = int(T4[cur])
            L_cyc = len(cyc)
            cyc_arr = np.array(cyc)
            cyc_rhs = rhs_scaled[cyc_arr]
            cyc_rhs -= cyc_rhs.sum() / L_cyc  # projecteer orthogonaal op ones
            V1_cyc = np.cumsum(cyc_rhs); V1_cyc -= V1_cyc.mean()
            V[cyc_arr] = V1_cyc
        V -= V.mean()
        return V

    rhs_min = rho_1 * ones - L1_min(ones)
    rhs_avg = rho_1 * ones - L1_avg(ones)

    V1_min = solve_via_cycle(rhs_min)
    V1_avg = solve_via_cycle(rhs_avg)

    std_min = float(np.std(V1_min))
    std_avg = float(np.std(V1_avg))

    corr_min_fd = float(np.corrcoef(V1_min, V1_fd)[0, 1]) if N > 10 else float('nan')
    corr_avg_fd = float(np.corrcoef(V1_avg, V1_fd)[0, 1]) if N > 10 else float('nan')

    print(f"  k={k}:")
    print(f"    Finite-diff V_1 (eps=0.005): std={std_fd:.5f}")
    print(f"    Min-linearisatie V_1:         std={std_min:.5f}  "
          f"corr(V1_min, fd)={corr_min_fd:.4f}")
    print(f"    Gem-linearisatie V_1:         std={std_avg:.5f}  "
          f"corr(V1_avg, fd)={corr_avg_fd:.4f}")
    print(f"    Verhouding std_fd/std_min = {std_fd/std_min:.4f}  "
          f"(analytisch bewijs vereist = 0.789/0.934 = 0.844 als puur lineair)")

    return std_fd, std_min, std_avg


print(f"221: Gecorrigeerde perturbatietheorie  (lam={LAM})")
print("=" * 65)
print()
print("(A) Numerieke afgeleiden van CV(eps) bij eps->0:")
for k in (8, 10, 12):
    numerical_derivative(k)

print()
print("(B) Gecorrigeerde V_1 (min vs gem linearisatie):")
for k in (8, 10, 12):
    correct_V1(k)
    print()

print("done")
