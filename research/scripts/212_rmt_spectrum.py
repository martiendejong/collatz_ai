"""
212_rmt_spectrum.py
===================
Methode 3: Willekeurige Matrix Theorie universaliteitsklasse van M = P_W L P_W.

Voor kleine k (k=8..11) berekenen we het volledige spectrum van M.
Dan vergelijken we de eigenwaardedichtheid met:
  - Wigner halve cirkel (GOE/GUE universaliteitsklasse)
  - Marchenko-Pastur verdeling (Wishart-klasse)
  - Uniforme/Poisson spacing (geen universaliteit)

Level spacing statistiek:
  - s_i = (lambda_{i+1} - lambda_i) / gemiddelde_spacing
  - Wigner-Dyson: P(s) ~ s * exp(-pi*s^2/4)   [GOE]
  - Poisson:      P(s) ~ exp(-s)               [geen correlaties]

Als M in de Wigner-Dyson klasse valt, is er een diepere willekeurige
matrix structuur achter de K-L-operator.

De dead-flat constante sigma_W/rho = 0.755 zou dan een RMT-universeel
getal kunnen zijn (vergelijkbaar met de GOE bulk-eigenwaarde-spreiding).
"""
import numpy as np
from math import log2
from scipy import linalg

ALPHA = log2(3.0)
LAM   = 1.70
A     = LAM ** -2.0
B1    = LAM ** (ALPHA - 2.0)
B3    = LAM ** (ALPHA - 1.0)


def build_M_matrix(k):
    """Bouw M = P_W o L o P_W expliciet als dense matrix voor kleine k."""
    N  = 3 ** (k - 1)
    i  = np.arange(N, dtype=np.int64)
    T4 = (4 * i + 2) % N
    s, r = np.divmod(i, 3)
    Nl = N // 3
    R1 = (4 * s) % Nl
    R3 = (2 * s + 1) % Nl
    m0, m2 = (r == 0), (r == 2)

    # Perron-vector voor frozen argmins
    v = np.ones(N, dtype=np.float64)
    for _ in range(300):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        w  = A * v[T4]
        w[m2] += B3 * cb[R3[m2]]
        w[m0] += B1 * cb[R1[m0]]
        v = w / w.max()
    # Frozen argmins
    stack = np.stack([v[:Nl], v[Nl:2*Nl], v[2*Nl:]])
    sel   = stack.argmin(axis=0).astype(np.int64)
    tgt2  = R3 + sel[R3] * Nl
    tgt0  = R1 + sel[R1] * Nl

    def PW(x):
        """Project op nulmiddel per top-digit drietal."""
        m = (x[:Nl] + x[Nl:2*Nl] + x[2*Nl:]) / 3.0
        y = x.copy()
        y[:Nl] -= m; y[Nl:2*Nl] -= m; y[2*Nl:] -= m
        return y

    def L_op(d):
        """Gelinieariseerde K-L-operator aan het vaste punt."""
        y = A * d[T4]
        y[m2] += B3 * d[tgt2[m2]]
        y[m0] += B1 * d[tgt0[m0]]
        return y

    # Bouw M als expliciete matrix: M e_j = P_W(L(P_W(e_j)))
    M = np.zeros((N, N), dtype=np.float64)
    e = np.zeros(N, dtype=np.float64)
    for j in range(N):
        e[j] = 1.0
        M[:, j] = PW(L_op(PW(e)))
        e[j] = 0.0

    return M


def level_spacing(eigenvalues):
    """Genormaliseerde level spacing s_i voor reele eigenwaarden."""
    ev = np.sort(np.real(eigenvalues))
    spacings = np.diff(ev)
    mean_s = np.mean(spacings)
    if mean_s < 1e-12:
        return np.array([])
    return spacings / mean_s


def wigner_dyson_p(s):
    """GOE Wigner-Dyson verdeling: P(s) = (pi/2)*s*exp(-pi*s^2/4)."""
    return (np.pi / 2) * s * np.exp(-np.pi * s**2 / 4)


def poisson_p(s):
    """Poisson: P(s) = exp(-s)."""
    return np.exp(-s)


print(f"Methode 3: RMT-spectrum van M = P_W L P_W  (lam={LAM})", flush=True)
print("=" * 65, flush=True)

for k in (6, 7, 8, 9):
    N = 3 ** (k - 1)
    print(f"\n--- k={k}  (N={N}) ---", flush=True)

    M = build_M_matrix(k)

    # Volledig spectrum
    eigs = linalg.eigvals(M)
    eigs_real = np.real(eigs)
    eigs_imag = np.imag(eigs)
    frac_real = float(np.mean(np.abs(eigs_imag) < 1e-8 * np.abs(eigs_real + 1e-12)))

    print(f"  Matrixgrootte: {N}x{N}", flush=True)
    print(f"  Spectrum: max_re={eigs_real.max():.5f}  "
          f"min_re={eigs_real.min():.5f}", flush=True)
    print(f"  Fractie reele eigenwaarden: {frac_real:.4f}", flush=True)

    # Level spacing van reele eigenwaarden
    re_eigs = eigs_real[np.abs(eigs_imag) < 0.01]
    if len(re_eigs) > 10:
        s = level_spacing(re_eigs)
        if len(s) > 0:
            # Vergelijk met GOE en Poisson via KS-statistiek
            from scipy.stats import ks_1samp, expon
            # KS test vs Poisson
            ks_poisson = ks_1samp(s, lambda x: 1 - np.exp(-x))
            print(f"  Level spacing: gem={s.mean():.4f}  "
                  f"std={s.std():.4f}  "
                  f"KS-vs-Poisson={ks_poisson.statistic:.4f} "
                  f"(p={ks_poisson.pvalue:.4f})", flush=True)

    # Eigenwaardedichtheid: histogram over [-sigma_W, sigma_W]
    sw = float(np.std(eigs_real))
    print(f"  std(Re eigenwaarden) = {sw:.5f}", flush=True)
    print(f"  sigma_W (power method) / vergelijking:", flush=True)

    # Top-3 eigenwaarden (buiten de nul-mode)
    abs_eigs = np.abs(eigs)
    top3 = np.sort(abs_eigs)[::-1][:5]
    print(f"  Top-5 |eigenwaarden|: {top3.round(5)}", flush=True)

    # Bulkstatistiek: Wigner halve cirkel fit
    bulk = eigs_real[np.abs(eigs_real) < 2*sw]
    bins = np.linspace(-2*sw, 2*sw, 20)
    hist, _ = np.histogram(bulk, bins=bins, density=True)
    # Wigner halve cirkel: rho(x) = (2/pi*R^2)*sqrt(R^2-x^2)
    R = sw * np.sqrt(2)
    x_mid = (bins[:-1] + bins[1:]) / 2
    wigner = (2 / (np.pi * R**2)) * np.sqrt(np.maximum(R**2 - x_mid**2, 0))
    residual = np.mean((hist - wigner)**2)
    print(f"  Wigner halve-cirkel fit MSE: {residual:.6f}", flush=True)

print("\ndone", flush=True)
