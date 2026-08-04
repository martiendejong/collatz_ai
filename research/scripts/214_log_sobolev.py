"""
214_log_sobolev.py
==================
Methode 5: Log-Sobolev-ongelijkheid voor de K-L-operator.

Een log-Sobolev-ongelijkheid (LSI) voor de operator L:
  Ent_mu(f) := E_mu[f log f] - E_mu[f] log E_mu[f]
             <= (1/2) * rho_LS * E(f, f)

waarbij E(f, f) = <(I - L/rho) f, f>_mu de Dirichlet-vorm is.

De LSI-constante rho_LS > 0 impliceert:
  - Hypercontractiviteit van e^{tL} (warmtekernsamenk)
  - Exponentiële convergentie naar de stationnairen toestand
  - Directe kwantitatieve ondergrens op CV

Numerieke aanpak:
  - Meet rho_LS via de variatieprincipe:
    rho_LS = min_f { Ent_mu(f^2) / E(f^2, f^2) }  [Bakry-Emery definitie]
  - Schat via Monte Carlo op willekeurige testfuncties f
  - Vergelijk met Poincaré constante rho_P = spectraalkloof (= 1 - sigma_W/rho)

Als rho_LS >= c * rho_P voor een universele c, dan is de operator
"hypercontractief" en geeft dat de sterkste versie van de CV-ondergrens.

Maat mu = de Perron-maat (genormaliseerde eigenvector).
"""
import numpy as np
from math import log2

ALPHA = log2(3.0)
LAM   = 1.70
A     = LAM ** -2.0
B1    = LAM ** (ALPHA - 2.0)
B3    = LAM ** (ALPHA - 1.0)


def build_system(k, n_iter=300):
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
    cb   = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
    w2   = A * v[T4]
    w2[m2] += B3 * cb[R3[m2]]
    w2[m0] += B1 * cb[R1[m0]]
    rho  = float(w2.sum() / v.sum())
    v   /= v.mean()
    mu   = v / v.sum()

    # Frozen argmins voor gelinieariseerde operator
    stack = np.stack([v[:Nl], v[Nl:2*Nl], v[2*Nl:]])
    sel   = stack.argmin(axis=0).astype(np.int64)
    tgt2  = R3 + sel[R3] * Nl
    tgt0  = R1 + sel[R1] * Nl

    def L_op(f):
        """Genormaliseerde gelinieariseerde operator (1/rho) * L * f."""
        y = A * f[T4]
        y[m2] += B3 * f[tgt2[m2]]
        y[m0] += B1 * f[tgt0[m0]]
        return y / rho

    return mu, v, rho, L_op, N


def entropy(f, mu):
    """Ent_mu(f) = E[f log f] - E[f] * log(E[f]) (f >= 0)."""
    Ef  = float(np.dot(f, mu))
    Efl = float(np.dot(f * np.log(f + 1e-300), mu))
    return Efl - Ef * np.log(Ef + 1e-300)


def dirichlet(f, mu, L_op):
    """
    Dirichlet-vorm E(f, f) = <f, (I - L) f>_mu
    = E[f^2] - E[f * L(f)].
    """
    Lf  = L_op(f)
    Ef2 = float(np.dot(f * f, mu))
    EfLf = float(np.dot(f * Lf, mu))
    return Ef2 - EfLf


def poincare_ratio(f, mu, L_op):
    """
    Poincaré ratio: Var_mu(f) / E(f, f).
    Poincaré constante = min over niet-constante f.
    """
    Ef  = float(np.dot(f, mu))
    Var = float(np.dot((f - Ef)**2, mu))
    D   = dirichlet(f - Ef, mu, L_op)
    if D < 1e-12:
        return np.inf
    return Var / D


def lsi_ratio(f, mu, L_op):
    """
    LSI ratio: Ent_mu(f^2) / E(f^2, f^2).
    LSI constante = min over niet-constante f.
    """
    f2  = f * f
    Ent = entropy(f2, mu)
    D   = dirichlet(f2, mu, L_op)
    if D < 1e-12 or Ent < 1e-12:
        return np.inf
    return Ent / D


print(f"Methode 5: Log-Sobolev-ongelijkheid  (lam={LAM})", flush=True)
print("=" * 65, flush=True)

n_trials = 500
rng = np.random.default_rng(42)

for k in (12, 13, 14):
    print(f"\n--- k={k}  (N={3**(k-1):,}) ---", flush=True)
    mu, v, rho, L_op, N = build_system(k)

    # Bereken Poincaré-constante via spectrale kloof (bekende waarde)
    # rho_P = 1 - sigma_W/rho (uit Obs 405: sigma_W/rho = 0.755)
    rho_P_theory = 1.0 - 0.755  # = 0.245

    poincare_vals = []
    lsi_vals = []

    for _ in range(n_trials):
        # Willekeurige testfunctie: lineaire combinatie van Perron-vector
        # en willekeurige componenten
        f_raw = rng.standard_normal(N)
        # Normeer f positief (voor LSI: f^2 > 0 altijd, maar we testen f direct)
        f = np.exp(0.3 * f_raw)  # altijd positief, niet te extreem
        f = f / np.dot(f, mu)    # genormaliseerd: E_mu[f] = 1

        p = poincare_ratio(f, mu, L_op)
        l = lsi_ratio(f, mu, L_op)
        if np.isfinite(p):
            poincare_vals.append(p)
        if np.isfinite(l):
            lsi_vals.append(l)

    if poincare_vals:
        rho_P_meas = float(np.min(poincare_vals))
        print(f"  Poincaré constante (gemeten): {rho_P_meas:.5f}", flush=True)
        print(f"  Poincaré constante (theorie 1-sigma_W/rho): {rho_P_theory:.5f}",
              flush=True)

    if lsi_vals:
        rho_LS = float(np.min(lsi_vals))
        print(f"  Log-Sobolev constante (gemeten): {rho_LS:.5f}", flush=True)
        if poincare_vals:
            ratio = rho_LS / rho_P_meas if rho_P_meas > 0 else float('nan')
            print(f"  rho_LS / rho_P = {ratio:.4f}  "
                  f"(>1 = log-Sobolev sterker dan Poincaré)", flush=True)

    # Implicatie voor CV: als LSI geldt met constante C,
    # dan Ent_mu(v) <= (1/2C) * E(v,v)
    # Via Bakry-Emery: implicatie voor eigenvector spreiding
    Ent_v  = entropy(v, mu)
    D_v    = dirichlet(v, mu, L_op)
    if D_v > 1e-12:
        implied_rho_LS = Ent_v / D_v if Ent_v > 0 else float('inf')
        print(f"  Ent_mu(v) = {Ent_v:.5f}   E(v,v) = {D_v:.5f}   "
              f"ratio = {implied_rho_LS:.4f}", flush=True)

print("\ndone", flush=True)
