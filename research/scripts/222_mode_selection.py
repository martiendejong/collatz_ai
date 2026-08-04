"""
222_mode_selection.py
=====================
Fourier-modeselectie: welke eigenmode van L_0 = A*T4 wordt geselecteerd
door de liftstoring L_1 bij kleine eps?

Bij eps=0: L_0 heeft N eigenvectoren u_j(i) = exp(2*pi*i*j*pos(i)/N)
waarbij pos(i) de positie van i in de T4-cyclus is.
Alle N eigenwaarden zijn A*exp(2*pi*i*j/N).

De liftstoring L_1 selecteert de mode j* die de hoogste
"Rayleigh-quotient versterking" heeft:
  R(j) = Re(< u_j, L_1 u_j >) / ||u_j||^2

De mode j* met maximale R(j) bepaalt de richting van de
Perron-eigenvector bij kleine eps.

Doelstelling:
  (A) Bereken R(j) voor alle j (of de dominante j's)
  (B) Identificeer j* = argmax R(j)
  (C) Verifieer dat de verwachte richting u_{j*} correleert met V1_fd
  (D) Bereken de theoretische CV(eps)/eps = ||u_{j*}|| / N voor kleine eps
"""
import numpy as np
from math import log2

ALPHA = log2(3.0)
LAM   = 1.70
A_coef  = LAM ** -2.0
B1_coef = LAM ** (ALPHA - 2.0)
B3_coef = LAM ** (ALPHA - 1.0)


def build_structure(k):
    N  = 3 ** (k - 1)
    i  = np.arange(N, dtype=np.int64)
    T4 = (4 * i + 2) % N
    s, r = np.divmod(i, 3)
    Nl = N // 3
    m0, m2 = (r == 0), (r == 2)
    R1 = (4 * s) % Nl
    R3 = (2 * s + 1) % Nl
    return T4, R1, R3, m0, m2, Nl, N


def get_cycle_positions(T4, N):
    """Vind de positie van elk element in de T4-cyclus."""
    pos = np.zeros(N, dtype=np.int64)
    cur = 0
    for p in range(N):
        pos[cur] = p
        cur = int(T4[cur])
    return pos


def L1_action(u, R1, R3, m0, m2, Nl):
    """L_1 u = B3*min(...) voor r=2, B1*min(...) voor r=0."""
    N = u.size
    cb = np.minimum(np.minimum(u[:Nl], u[Nl:2*Nl]), u[2*Nl:])
    w  = np.zeros_like(u)
    w[m2] = B3_coef * cb[R3[m2]]
    w[m0] = B1_coef * cb[R1[m0]]
    return w


def rayleigh_quotient(j, pos, R1, R3, m0, m2, Nl, N):
    """
    Rayleigh quotient R(j) = Re(< u_j, L_1 u_j >) voor mode j.
    u_j(i) = exp(2*pi*i*j*pos(i)/N)
    """
    theta = 2.0 * np.pi * j * pos / N
    uj = np.cos(theta) + 1j * np.sin(theta)
    # L_1 u_j (complexe versie): min van complexe getallen = min van Re
    # Correctie: de echte operatorverking gebruikt min over ECHTE delen
    # Maar u_j is complex; we werken met Re(u_j) als echte perturbatie
    uj_re = np.cos(theta).copy()
    L1_uj = L1_action(uj_re, R1, R3, m0, m2, Nl)
    # Rayleigh quotient voor echte eigenmode:
    R = float(np.dot(uj_re, L1_uj)) / float(np.dot(uj_re, uj_re))
    return R


def find_dominant_mode(k, max_modes=100):
    """Vind de dominante Fourier-mode van L_1 op de eigenruimte van L_0."""
    T4, R1, R3, m0, m2, Nl, N = build_structure(k)
    pos = get_cycle_positions(T4, N)

    print(f"  k={k}  N={N}")

    # Bereken Rayleigh-quotienten voor alle modes j=0..max_modes
    max_j = min(max_modes, N // 2)
    Rs = []
    for j in range(max_j + 1):
        R = rayleigh_quotient(j, pos, R1, R3, m0, m2, Nl, N)
        Rs.append((j, R))

    Rs.sort(key=lambda x: -x[1])
    print(f"  Top-10 Rayleigh-quotienten:")
    for j, R in Rs[:10]:
        print(f"    j={j:5d}:  R={R:.6f}")

    j_star = Rs[0][0]
    R_star = Rs[0][1]

    # Bereken de eigenvector voor j*
    theta_star = 2.0 * np.pi * j_star * pos / N
    u_star = np.cos(theta_star)
    u_star -= u_star.mean()

    # Vergelijk met de numerieke V1_fd (perron-eigenvector bij kleine eps)
    # (Berekening van v(eps) voor kleine eps)
    eps_test = 0.005
    v = np.ones(N, dtype=np.float64)
    B1 = eps_test * B1_coef
    B3 = eps_test * B3_coef
    for _ in range(500):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        w  = A_coef * v[T4]
        w[m2] += B3 * cb[R3[m2]]
        w[m0] += B1 * cb[R1[m0]]
        v = w / w.max()
    v /= v.mean()
    V1_fd = v - 1.0
    V1_fd -= V1_fd.mean()

    corr_star = float(np.corrcoef(u_star, V1_fd)[0, 1]) if np.std(u_star) > 1e-10 else float('nan')
    print(f"  Mode j*={j_star}: R*={R_star:.6f}  corr(u_j*, V1_fd)={corr_star:.4f}")

    # Bereken ook for j=1 en j=N-1 (laagste niet-triviale modes)
    R_1 = rayleigh_quotient(1, pos, R1, R3, m0, m2, Nl, N)
    R_N1 = rayleigh_quotient(N-1, pos, R1, R3, m0, m2, Nl, N)
    print(f"  Mode j=1:    R={R_1:.6f}")
    print(f"  Mode j=N-1:  R={R_N1:.6f}")
    print(f"  Mode j=0 (uniform): R={rayleigh_quotient(0, pos, R1, R3, m0, m2, Nl, N):.6f}")

    # Theoretisch: als V1_fd = c*u_{j*}, dan CV(eps)/eps ~ std(u_{j*})
    std_u = float(np.std(u_star))
    print(f"  std(u_{j_star}) = {std_u:.5f}  (theoretisch: CV(eps)/eps bij kleine eps)")
    print()
    return j_star, R_star, corr_star


print(f"222: Fourier-modeselectie  (lam={LAM})")
print("=" * 65)

for k in (5, 6, 7, 8):
    find_dominant_mode(k)

print("done")
