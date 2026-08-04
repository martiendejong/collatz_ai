"""
223_deadflat_catalog.py
=======================
Catalogus van alle dead-flat constanten van de K-L-operator (lam=1.70).

Een dead-flat constante is k-invariant in de limiet k->inf.
Bewijs-relevantie: als de waarde positief (of < 1) is, is de bijbehorende
bewijs-stap k-uniform gesloten.
"""
import numpy as np
from math import log2

ALPHA = log2(3.0)
LAM   = 1.70
A  = LAM ** -2.0
B1 = LAM ** (ALPHA - 2.0)
B3 = LAM ** (ALPHA - 1.0)


def compute_all(k, n_iter=300):
    N  = 3 ** (k - 1)
    i  = np.arange(N, dtype=np.int64)
    T4 = (4 * i + 2) % N
    s, r = np.divmod(i, 3)
    Nl = N // 3
    m0, m2 = (r == 0), (r == 2)
    R1 = (4 * s) % Nl
    R3 = (2 * s + 1) % Nl

    # Perron-eigenvector (genormaliseerd op gemiddelde = 1)
    v = np.ones(N, dtype=np.float64)
    for _ in range(n_iter):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        w  = A * v[T4]
        w[m2] += B3 * cb[R3[m2]]
        w[m0] += B1 * cb[R1[m0]]
        v = w / w.max()
    v /= v.mean()   # normaliseer op gemiddelde = 1

    # Eigenwaarde rho
    cb2 = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
    w2  = A * v[T4]
    w2[m2] += B3 * cb2[R3[m2]]
    w2[m0] += B1 * cb2[R1[m0]]
    rho = float(w2.mean())   # w2.mean() = rho * v.mean() = rho * 1 = rho

    # CV
    cv = float(np.std(v))    # std/mean = std (want mean=1)

    # sigma_W / rho via power-iteratie op P_W L P_W
    stack = np.stack([v[:Nl], v[Nl:2*Nl], v[2*Nl:]])
    sel   = stack.argmin(axis=0).astype(np.int64)
    tgt2  = R3 + sel[R3] * Nl
    tgt0  = R1 + sel[R1] * Nl

    def PW(x):
        m = (x[:Nl] + x[Nl:2*Nl] + x[2*Nl:]) / 3.0
        y = x.copy(); y[:Nl] -= m; y[Nl:2*Nl] -= m; y[2*Nl:] -= m
        return y

    rng = np.random.default_rng(42)
    d   = PW(rng.standard_normal(N)); d /= np.linalg.norm(d)
    rates = []
    for _ in range(200):
        y = A * d[T4]
        y[m2] += B3 * d[tgt2[m2]]
        y[m0] += B1 * d[tgt0[m0]]
        y = PW(y); nrm = np.linalg.norm(y); rates.append(nrm); d = y/(nrm+1e-300)
    sw_rho = float(np.exp(np.mean(np.log(np.array(rates[-80:])+1e-300)))) / rho

    # Herschakelpercentage (fractie van type-2 knopen waarbij de argmin wisselt)
    # sel[i] = 0: D1-component is minimum, 2: D2-component, etc.
    n_switch = int((sel[R3[m2]] > 0).sum())   # R3 verwijst naar laag-0 pos; >0 = andere sel
    switch_pct = 100.0 * n_switch / (m2.sum() if m2.sum() > 0 else 1)

    # v2/v1 ratio (type-2 gemiddeld t.o.v. type-1)
    v_type2 = float(v[m2].mean())
    v_type1 = float(v[m0].mean())   # type-0 = D1 in notatie
    v2_v1 = v_type2 / v_type1 if v_type1 > 0 else float('nan')

    # Log-log helling (helling van log(sorted v) vs log(sorted rank))
    # Vergelijkbaar met Obs 408 "log-log slope"
    v_sorted = np.sort(v)
    ranks = np.arange(1, N+1) / N
    ok = (v_sorted > 0.01) & (ranks > 0.05) & (ranks < 0.95)
    if ok.sum() >= 20:
        slope_ll, _ = np.polyfit(np.log(ranks[ok]), np.log(v_sorted[ok]), 1)
    else:
        slope_ll = float('nan')

    return {
        'k': k, 'N': N, 'rho': rho, 'cv': cv,
        'sw_rho': sw_rho, 'switch_pct': switch_pct,
        'v2_v1': v2_v1, 'slope_ll': slope_ll
    }


print(f"223: Catalogus dead-flat constanten  (lam={LAM})")
print("=" * 75)
print()

print("=== ANALYTISCH k-INVARIANTE CONSTANTEN ===")
print(f"  rho_1 = (B3+B1)/3 = {(B3+B1)/3:.6f}  [exact: lam^(al-1)/3 + lam^(al-2)/3]")
print()

print("=== NUMERIEK GEMETEN DEAD-FLAT CONSTANTEN ===")
print(f"{'k':>3}  {'CV':>7}  {'dCV':>7}  {'sw/rho':>7}  {'sw%':>6}  "
      f"{'v2/v1':>6}  {'LL-sl':>6}")

prev_cv = None
for k in range(10, 16):
    r = compute_all(k)
    dCV = (r['cv'] - prev_cv) if prev_cv is not None else float('nan')
    print(f"  {k:2d}  {r['cv']:7.5f}  {dCV:7.5f}  {r['sw_rho']:7.5f}  "
          f"{r['switch_pct']:5.2f}%  {r['v2_v1']:6.3f}  {r['slope_ll']:6.4f}")
    prev_cv = r['cv']

print()
print("=== COMPLETE CATALOGUS (10 dead-flat constanten) ===")
lines = [
    "  #  Grootheid                     Waarde      Obs     Bewijs-relevantie",
    "  " + "-" * 68,
    "  1  sw/rho (sigma_W/rho)          0.755      Obs405  spectrale kloof 0.245",
    "  2  r_real/rho                    0.491      Obs406  spectr. kloof herbevestigd",
    "  3  Herschakelpercentage          88.90%     Obs404  33/37 per stap",
    "  4  Log-log helling               1.2450     Obs408  koppelingssterkte stabiel",
    "  5  CV-increment per k            0.001/stap Obs413  CV convergeert",
    "  6  std(Re eigenw. M)             0.251      Obs416  nieuwe universaliteitsklasse",
    "  7  rho_LS (Log-Sobolev)          0.855      Obs419  hypercontractiviteit",
    "  8  Gesorteerde W1 (Wasserstein)  0.273      Obs418  Perron-maat convergentie",
    "  9  rho_1 = (B3+B1)/3            0.7221     Obs425  analytisch exact (k-vrij)",
    "  10 CV(eps)/eps bij eps=0.01      0.9225     Obs426  zelf-similariteit eigenvec.",
]
print("\n".join(lines))
print("done")
