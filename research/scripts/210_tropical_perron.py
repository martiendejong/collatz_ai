"""
210_tropical_perron.py
======================
Methode 1: Tropische algebra van de min-smoothing operator.

De K-L-iteratie bevat cb = min(v1,v2,v3) -- een tropische operatie.
In log-ruimte w_i = log v_i wordt de operator:

  w_i <- log( A*exp(w_{T4}) + B3*exp(min_j w_{R3+j*Nl}) )   [type 2]
  w_i <- log( A*exp(w_{T4}) )                                 [type 1]
  w_i <- log( A*exp(w_{T4}) + B1*exp(min_j w_{R1+j*Nl}) )   [type 0]

De TROPISCHE benadering (T->0 limiet, max-min algebra):
  w_i <- max( log A + w_{T4},  log B3 + min_j w_{R3+j*Nl} ) [type 2]
  w_i <- log A + w_{T4}                                       [type 1]
  w_i <- max( log A + w_{T4},  log B1 + min_j w_{R1+j*Nl} ) [type 0]

Dit is een max-min operator. Zijn vaste punt = tropische Perron-vector.
De tropische eigenwaarde = maximaal gemiddeld cyclusgewicht in de digraph.

Vergelijking tropisch vs. echt:
- Als ze dicht bij elkaar liggen -> de dominante term bestuurt de eigenvector
- De kloof = kwantificering van de "endlichtetemperatuur"-correctie
- CV van de tropische vector = ondergrens op CV van de echte vector
"""
import numpy as np
from math import log2, log, exp

ALPHA = log2(3.0)
LAM   = 1.70
A     = LAM ** -2.0
B1    = LAM ** (ALPHA - 2.0)
B3    = LAM ** (ALPHA - 1.0)
logA  = log(A)
logB1 = log(B1)
logB3 = log(B3)


def build(k):
    N  = 3 ** (k - 1)
    i  = np.arange(N, dtype=np.int64)
    T4 = (4 * i + 2) % N
    s, r = np.divmod(i, 3)
    Nl = N // 3
    R1 = (4 * s) % Nl
    R3 = (2 * s + 1) % Nl
    m0, m1, m2 = (r == 0), (r == 1), (r == 2)
    return N, Nl, T4, R1, R3, m0, m1, m2


def tropical_perron(k, n_iter=600):
    """Max-min iteratie in log-ruimte. Eigenwaarde via gemiddelde verschuiving."""
    N, Nl, T4, R1, R3, m0, m1, m2 = build(k)
    w = np.zeros(N, dtype=np.float64)  # log v, start uniform

    for it in range(n_iter):
        cb_log = np.minimum(
            np.minimum(w[:Nl], w[Nl:2*Nl]), w[2*Nl:]
        )  # min in log-ruimte = log(min v)

        w_new = np.empty(N, dtype=np.float64)
        # type 1: enkel walk-on-4
        w_new[m1] = logA + w[T4[m1]]
        # type 2: max van walk en B3-lift
        w_new[m2] = np.maximum(logA + w[T4[m2]],
                                logB3 + cb_log[R3[m2]])
        # type 0: max van walk en B1-lift
        w_new[m0] = np.maximum(logA + w[T4[m0]],
                                logB1 + cb_log[R1[m0]])

        # verschuif met gemiddelde = verwijder de eigenwaarde-component
        shift = w_new.mean()
        w = w_new - shift

    # tropische eigenwaarde = gemiddelde verschuiving per stap (laatste 100)
    shifts = []
    for _ in range(100):
        cb_log = np.minimum(np.minimum(w[:Nl], w[Nl:2*Nl]), w[2*Nl:])
        w_new = np.empty(N, dtype=np.float64)
        w_new[m1] = logA + w[T4[m1]]
        w_new[m2] = np.maximum(logA + w[T4[m2]], logB3 + cb_log[R3[m2]])
        w_new[m0] = np.maximum(logA + w[T4[m0]], logB1 + cb_log[R1[m0]])
        shift = w_new.mean()
        shifts.append(shift)
        w = w_new - shift

    trop_rho = exp(np.mean(shifts))
    return w, trop_rho


def real_perron_log(k, n_iter=300):
    """Gewone Perron-iteratie, teruggegeven in log-ruimte."""
    N, Nl, T4, R1, R3, m0, m1, m2 = build(k)
    v = np.ones(N, dtype=np.float64)
    for _ in range(n_iter):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        w2 = A * v[T4]
        w2[m2] += B3 * cb[R3[m2]]
        w2[m0] += B1 * cb[R1[m0]]
        v = w2 / w2.max()
    cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
    w2 = A * v[T4]
    w2[m2] += B3 * cb[R3[m2]]
    w2[m0] += B1 * cb[R1[m0]]
    rho = float(w2.sum() / v.sum())
    v /= v.mean()
    return np.log(v), rho


print(f"Methode 1: Tropische Perron-vector  (lam={LAM})", flush=True)
print("=" * 65, flush=True)

for k in (12, 13, 14, 15):
    print(f"\n--- k={k}  (N={3**(k-1):,}) ---", flush=True)

    w_trop, trop_rho = tropical_perron(k)
    w_real, real_rho = real_perron_log(k)

    # Beide centreren
    w_trop -= w_trop.mean()
    w_real -= w_real.mean()

    # CV in log-ruimte (std van log v)
    cv_trop = float(np.std(w_trop))
    cv_real = float(np.std(w_real))

    # Correlatie tropisch vs. echt
    corr = float(np.corrcoef(w_trop, w_real)[0, 1])

    # L2-afstand genormaliseerd
    diff = w_trop - w_real
    l2 = float(np.sqrt(np.mean(diff**2)))

    # Fractie nodes waar lift wint over walk (in tropische versie)
    N, Nl, T4, R1, R3, m0, m1, m2 = build(k)
    cb_log = np.minimum(np.minimum(w_trop[:Nl], w_trop[Nl:2*Nl]), w_trop[2*Nl:])
    walk2  = logA + w_trop[T4[m2]]
    lift2  = logB3 + cb_log[R3[m2]]
    frac_lift_wins = float((lift2 > walk2).mean())

    print(f"  rho_trop  = {trop_rho:.6f}   rho_real = {real_rho:.6f}",
          flush=True)
    print(f"  CV_trop (log-std) = {cv_trop:.5f}   CV_real = {cv_real:.5f}",
          flush=True)
    print(f"  corr(trop, real)  = {corr:.5f}", flush=True)
    print(f"  L2-verschil       = {l2:.5f}", flush=True)
    print(f"  lift wint (type2) = {frac_lift_wins:.4f} van alle type-2 knopen",
          flush=True)

    # Kwantielen van tropische vector
    q10 = float(np.quantile(np.exp(w_trop), 0.10))
    q90 = float(np.quantile(np.exp(w_trop), 0.90))
    print(f"  trop Q10/Q90      = {q10:.4f} / {q90:.4f}  (ratio {q90/q10:.2f}x)",
          flush=True)
    q10r = float(np.quantile(np.exp(w_real), 0.10))
    q90r = float(np.quantile(np.exp(w_real), 0.90))
    print(f"  real Q10/Q90      = {q10r:.4f} / {q90r:.4f}  (ratio {q90r/q10r:.2f}x)",
          flush=True)

print("\ndone", flush=True)
