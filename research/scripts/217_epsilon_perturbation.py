"""
217_epsilon_perturbation.py
===========================
Methode 8: Gestuurde perturbatie van de "+1"-correctie.

Vervang de Collatz-operator (3n+1)/2 door (3n+Îµ)/2 voor Îµ âˆˆ [0,1].

In de K-L-context: de "+1" verschijnt als de B3/B1-gewichten
(ze zijn Î»^{Î±-1} en Î»^{Î±-2} door de correctie-accumulatie).
De Îµ-versie correspondeert met gewichten:
  B3(Îµ) = Îµ * B3_orig + (1-Îµ) * B3_tropical
  B1(Îµ) = Îµ * B1_orig + (1-Îµ) * B1_tropical

waarbij B_tropical = 1 (de puur geometrische Bound Îµ=0).

Bij Îµ=0: puur geometrisch systeem (3n/2 in plaats van (3n+1)/2).
  Dan is n*3/2 nooit een geheel getal voor oneven n -> geen cycli
  analytisch bewijsbaar! Dit is de "gemakkelijke" Bound.
Bij Îµ=1: originele Collatz.

We meten hoe de dead-flat constanten en de koppeling veranderen als
Îµ varieert van 0 naar 1:
  - sigma_W(Îµ)/rho(Îµ)  (operator spreiding)
  - CV(v^(k), Îµ)       (eigenvector spreiding)
  - corr(logL, logS, Îµ) (rijkheid-ruwheid koppeling)
  - f2(Îµ)              (min-smoothing factor)

Als alles continu is in Îµ, geeft de perturbatieroute:
bewijs voor Îµ=0 + continuÃ¯teitsargument = bewijs voor Îµ=1.
"""
import numpy as np
from math import log2

ALPHA = log2(3.0)
LAM   = 1.70


def make_operator(lam, eps):
    """
    K-L operator met Îµ-perturbatie.
    B3(Îµ) = lam^{Î±-1+Îµ*Î´} waarbij Î´=0 het origineel is.

    Directere aanpak: interpoleer gewichten lineair.
    B_geo = 1.0  (zuiver geometrisch, geen correctie-accumulatie)
    B_col = lam^{Î±-1}  (originele Collatz)
    B3(Îµ) = (1-Îµ)*B_geo + Îµ*B_col

    Maar dit klopt niet precies: B3 = lam^{Î±-1} INCLUSIEF de +1 via
    de definitie van lam. Bij Îµ=0 willen we de "puur multiplicatieve"
    operator zonder additieve correctie.

    Praktische benadering: schaal de liftgewichten lineair.
    """
    A  = lam ** -2.0
    B3_full = lam ** (ALPHA - 1.0)
    B1_full = lam ** (ALPHA - 2.0)
    # Bij Îµ=0: enkel de walk-op-4 (geen lift) = triviale operator
    # Bij Îµ=1: volledige K-L
    B3 = eps * B3_full
    B1 = eps * B1_full
    return A, B1, B3


def perron_and_stats(k, lam, eps, n_iter=300, G=8):
    """Perron-vector en statistieken voor gegeven Îµ."""
    N  = 3 ** (k - 1)
    i  = np.arange(N, dtype=np.int64)
    T4 = (4 * i + 2) % N
    s, r = np.divmod(i, 3)
    Nl = N // 3
    m0, m2 = (r == 0), (r == 2)
    R1 = (4 * s) % Nl
    R3 = (2 * s + 1) % Nl
    A, B1, B3 = make_operator(lam, eps)

    v = np.ones(N, dtype=np.float64)
    for _ in range(n_iter):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        w  = A * v[T4]
        if eps > 0:
            w[m2] += B3 * cb[R3[m2]]
            w[m0] += B1 * cb[R1[m0]]
        v = w / w.max()
    cb2  = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
    w2   = A * v[T4]
    if eps > 0:
        w2[m2] += B3 * cb2[R3[m2]]
        w2[m0] += B1 * cb2[R1[m0]]
    rho  = float(w2.sum() / v.sum())
    v   /= v.mean()
    cv   = float(np.std(v) / np.mean(v))

    # sigma_W via power method op P_W L P_W
    stack = np.stack([v[:Nl], v[Nl:2*Nl], v[2*Nl:]])
    sel   = stack.argmin(axis=0).astype(np.int64)
    tgt2  = R3 + sel[R3] * Nl
    tgt0  = R1 + sel[R1] * Nl

    def PW(x):
        m = (x[:Nl] + x[Nl:2*Nl] + x[2*Nl:]) / 3.0
        y = x.copy()
        y[:Nl] -= m; y[Nl:2*Nl] -= m; y[2*Nl:] -= m
        return y

    rng = np.random.default_rng(1)
    d   = PW(rng.standard_normal(N))
    d  /= np.linalg.norm(d)
    rates = []
    for _ in range(200):
        y = A * d[T4]
        if eps > 0:
            y[m2] += B3 * d[tgt2[m2]]
            y[m0] += B1 * d[tgt0[m0]]
        y = PW(y)
        nrm = np.linalg.norm(y)
        rates.append(nrm)
        d = y / (nrm + 1e-300)
    sw_rho = float(np.exp(np.mean(np.log(np.array(rates[-80:]) + 1e-300)))) / rho

    # Rijkheid-ruwheid koppeling (vereenvoudigd: corr(level, spread) per kolom-blok)
    if N >= 27:
        Nl3 = Nl // 3
        M_mat = np.stack([v[:Nl], v[Nl:2*Nl], v[2*Nl:]])
        level  = M_mat.mean(axis=0)[:Nl3]
        spread = M_mat.std(axis=0)[:Nl3]
        ok = spread > 1e-12
        if ok.sum() > 10:
            corr_ls = float(np.corrcoef(np.log(level[ok]+1e-12),
                                         np.log(spread[ok]+1e-12))[0,1])
        else:
            corr_ls = float('nan')
    else:
        corr_ls = float('nan')

    return rho, cv, sw_rho, corr_ls


print(f"Methode 8: eps-perturbatie  (lam={LAM}, k=14)", flush=True)
print("=" * 65, flush=True)
print(f"  eps    rho      CV(v)    sw/rho   corr(logL,logS)", flush=True)

k = 14
for eps in np.linspace(0.0, 1.0, 11):
    rho, cv, sw_rho, corr_ls = perron_and_stats(k, LAM, eps)
    corr_str = 'nan' if np.isnan(corr_ls) else f'{corr_ls:.5f}'
    print(f"  {eps:.2f}   {rho:.5f}  {cv:.5f}  {sw_rho:.5f}  {corr_str}",
          flush=True)

print(f"\n  eps=0: puur walk-op-4, geen liften.", flush=True)
print(f"  eps=1: volledige K-L-operator.", flush=True)
print(f"  Continuiteit -> perturbatieroute voor het bewijs.", flush=True)

# Extra: schat helling d(CV)/d(eps) bij eps=1
eps_hi = 0.98
eps_lo = 0.92
_, cv_hi, _, _ = perron_and_stats(k, LAM, eps_hi)
_, cv_lo, _, _ = perron_and_stats(k, LAM, eps_lo)
slope = (cv_hi - cv_lo) / (eps_hi - eps_lo)
print(f"\n  d(CV)/d(eps) bij eps~1: {slope:.4f}", flush=True)

print("\ndone", flush=True)

