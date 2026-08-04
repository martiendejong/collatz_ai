"""
224_type_averages.py
====================
Typegemiddelden van de Perron-eigenvector van de K-L-operator.

Definitie:
  a = <v>_{r=0}  (D1-knopen: i mod 3 = 0, kleine liftterm B1)
  b = <v>_{r=2}  (D3-knopen: i mod 3 = 2, grote liftterm B3)
  c = <v>_{r=1}  (D2-knopen: i mod 3 = 1, GEEN liftterm)

Normalisatie: gemiddelde over ALLE N knopen = 1, dus a+b+c = 3.

ANALYTISCHE VERWACHTING:
T4: i -> (4i+2) mod N heeft  (4i+2) mod 3 = (i-1) mod 3.
  Dus T4 mapt r=0 -> r=2, r=2 -> r=1, r=1 -> r=0.
Voor r=1 knopen: rho*v[i] = A*v[T4(i)] met T4(i) in r=0.
Gemiddeld: rho*c = A*a  =>  c/a = A/rho   (EXACT).
Dus c is volledig bepaald door a, A en rho.

Vraag: is b/a ook analytisch bepaald?
"""
import numpy as np
from math import log2

ALPHA = log2(3.0)
LAM   = 1.70
A  = LAM ** -2.0
B1 = LAM ** (ALPHA - 2.0)
B3 = LAM ** (ALPHA - 1.0)


def analyze_type_averages(k, n_iter=400):
    N  = 3 ** (k - 1)
    i  = np.arange(N, dtype=np.int64)
    T4 = (4 * i + 2) % N
    s, r = np.divmod(i, 3)
    Nl = N // 3
    m0, m2, m1 = (r == 0), (r == 2), (r == 1)
    R1 = (4 * s) % Nl
    R3 = (2 * s + 1) % Nl

    # Perron-eigenvector (gemiddelde = 1)
    v = np.ones(N, dtype=np.float64)
    for _ in range(n_iter):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        w  = A * v[T4]
        w[m2] += B3 * cb[R3[m2]]
        w[m0] += B1 * cb[R1[m0]]
        v = w / w.max()
    v /= v.mean()

    # Eigenwaarde rho
    cb2 = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
    w2  = A * v[T4]
    w2[m2] += B3 * cb2[R3[m2]]
    w2[m0] += B1 * cb2[R1[m0]]
    rho = float(w2.mean())

    # Typegemiddelden
    a = float(v[m0].mean())  # D1 (r=0)
    b = float(v[m2].mean())  # D3 (r=2)
    c = float(v[m1].mean())  # D2 (r=1, geen lift)

    # Verificatie c/a = A/rho
    c_a_measured  = c / a
    c_a_predicted = A / rho

    # Verificatie a+b+c = 3
    sum_abc = a + b + c

    # Verificatie via directe gemiddelde van de r=1 vergelijking:
    # rho * v[i] = A * v[T4(i)] voor r=1
    # Gemiddeld: rho * c =? A * <v[T4(i)]>_{r=1}
    avg_vT4_r1 = float(v[T4[m1]].mean())  # T4(r=1) -> r=0, dus dit = <v>_{r=0} = a
    lhs_r1 = rho * c
    rhs_r1 = A * avg_vT4_r1

    # Verificatie via directe gemiddelde van de r=0 vergelijking:
    # rho * v[i] = A * v[T4(i)] + B1 * cb(R1(s)) voor r=0
    avg_vT4_r0    = float(v[T4[m0]].mean())  # T4(r=0) -> r=2
    avg_cb_R1     = float(np.minimum(v[R1[m0]],
                                     np.minimum(v[R1[m0] + Nl],
                                                v[R1[m0] + 2*Nl])).mean())
    lhs_r0 = rho * a
    rhs_r0 = A * avg_vT4_r0 + B1 * avg_cb_R1

    # Verificatie via directe gemiddelde van de r=2 vergelijking:
    avg_vT4_r2    = float(v[T4[m2]].mean())  # T4(r=2) -> r=1
    avg_cb_R3     = float(np.minimum(v[R3[m2]],
                                     np.minimum(v[R3[m2] + Nl],
                                                v[R3[m2] + 2*Nl])).mean())
    lhs_r2 = rho * b
    rhs_r2 = A * avg_vT4_r2 + B3 * avg_cb_R3

    return {
        'k': k, 'N': N, 'rho': rho,
        'a': a, 'b': b, 'c': c,
        'b_a': b/a, 'c_a': c/a, 'c_a_pred': c_a_predicted,
        'sum_abc': sum_abc,
        'lhs_r1': lhs_r1, 'rhs_r1': rhs_r1,   # should match exactly (no lift)
        'lhs_r0': lhs_r0, 'rhs_r0': rhs_r0,
        'lhs_r2': lhs_r2, 'rhs_r2': rhs_r2,
        'avg_vT4_r0': avg_vT4_r0,   # = b (T4 maps r=0 -> r=2)
        'avg_vT4_r1': avg_vT4_r1,   # = a (T4 maps r=1 -> r=0)
        'avg_vT4_r2': avg_vT4_r2,   # = c (T4 maps r=2 -> r=1)
        'avg_cb_R1': avg_cb_R1,
        'avg_cb_R3': avg_cb_R3,
    }


print(f"224: Typegemiddelden K-L Perron-eigenvector  (lam={LAM})")
print("=" * 70)
print()
print(f"Analytische verwachting: c/a = A/rho = {A:.6f}/rho")
print()

print(f"{'k':>3}  {'rho':>7}  {'a=D1':>7}  {'b=D3':>7}  {'c=D2':>7}  "
      f"{'b/a':>6}  {'c/a_m':>7}  {'c/a_p':>7}  {'err':>9}")

results = []
for k in range(8, 16):
    r = analyze_type_averages(k)
    results.append(r)
    err_c_a = abs(r['c_a'] - r['c_a_pred'])
    print(f"  {k:2d}  {r['rho']:7.5f}  {r['a']:7.5f}  {r['b']:7.5f}  {r['c']:7.5f}  "
          f"{r['b_a']:6.4f}  {r['c_a']:7.5f}  {r['c_a_pred']:7.5f}  {err_c_a:.2e}")

print()
print("=== VERIFICATIE GEMIDDELDE VASTE-PUNT-VERGELIJKINGEN ===")
for r in results:
    err_r1 = abs(r['lhs_r1'] - r['rhs_r1'])
    err_r0 = abs(r['lhs_r0'] - r['rhs_r0'])
    err_r2 = abs(r['lhs_r2'] - r['rhs_r2'])
    print(f"  k={r['k']:2d}: r=1 err={err_r1:.2e}  r=0 err={err_r0:.2e}  r=2 err={err_r2:.2e}")

print()
print("=== T4-PERMUTATIE-VERIFICATIE ===")
for r in results[:4]:
    # T4 maps r=0->r=2, r=2->r=1, r=1->r=0
    # Dus avg_vT4_r0 should = b, avg_vT4_r1 should = a, avg_vT4_r2 should = c
    print(f"  k={r['k']:2d}: <vT4>_r0={r['avg_vT4_r0']:.5f} vs b={r['b']:.5f} "
          f"| <vT4>_r1={r['avg_vT4_r1']:.5f} vs a={r['a']:.5f} "
          f"| <vT4>_r2={r['avg_vT4_r2']:.5f} vs c={r['c']:.5f}")

print()
print("=== ANALYTISCHE SYSTEM VAN VERGELIJKINGEN (sluitende vgl voor b/a) ===")
print()
print("Gedefinieerd: a = <v>_D1, b = <v>_D3, c = <v>_D2, <cb_R1> = theta_1, <cb_R3> = theta_3")
print("Stelsel:")
print("  rho*c = A*a                    (exact, r=1 heeft geen lift)")
print("  rho*a = A*b + B1*theta_1       (gemiddeld over r=0)")
print("  rho*b = A*c + B3*theta_3       (gemiddeld over r=2)")
print()
print("Eliminatie: c = (A/rho)*a, dan:")
print("  rho*a = A*b + B1*theta_1  =>  b = (rho*a - B1*theta_1)/A")
print("  rho*b = A*(A/rho)*a + B3*theta_3  =>  rho*b = (A^2/rho)*a + B3*theta_3")
print()
print("Als b/a = const, dan theta_3/a = const ook. Dat wil zeggen:")
print("  <cb_R3> / <v>_D1 is k-invariant.")
print()
print("Gemeten theta_3/a en theta_1/a:")
for r in results:
    theta3_a = r['avg_cb_R3'] / r['a']
    theta1_a = r['avg_cb_R1'] / r['a']
    print(f"  k={r['k']:2d}: theta_3/a={theta3_a:.5f}  theta_1/a={theta1_a:.5f}")

print()
print("=== ANALYTISCHE UITDRUKKING VOOR b/a ===")
print()
# From: rho*a = A*b + B1*theta_1  =>  b = (rho*a - B1*theta_1)/A
# => b/a = (rho - B1*(theta_1/a)) / A
# From: rho*b = (A^2/rho)*a + B3*theta_3
# => rho*(b/a) = A^2/rho + B3*(theta_3/a)
# => b/a = A^2/rho^2 + (B3/rho)*(theta_3/a)
print("Beide uitdrukkingen voor b/a (consistentiecheck):")
for r in results:
    theta3_a = r['avg_cb_R3'] / r['a']
    theta1_a = r['avg_cb_R1'] / r['a']
    ba_via_r0 = (r['rho'] - B1 * theta1_a) / A   # uit r=0 vgl
    ba_via_r2 = A**2 / r['rho']**2 + (B3 / r['rho']) * theta3_a  # uit r=2 vgl
    print(f"  k={r['k']:2d}: via r=0: {ba_via_r0:.5f}  via r=2: {ba_via_r2:.5f}  "
          f"direct: {r['b_a']:.5f}")

print()
print("done")
