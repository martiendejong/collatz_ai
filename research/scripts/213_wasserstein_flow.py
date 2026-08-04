"""
213_wasserstein_flow.py
=======================
Methode 4: Wasserstein-gradiëntstroom van de Perron-maatreeks.

De Perron-vector v^(k) definieert een kansmaat mu_k op {0,...,N-1}.
De Wasserstein-afstand W_1(mu_k, mu_{k+1}) meet hoe ver de maat
beweegt van diepte k naar k+1.

Intuïtie: als mu_k een gradiëntstroom is van een functionaal F,
dan neemt F monotoon af. De snelheid van afname = |grad F|^2.

We meten:
  (A) W_1(mu_k, mu_{k+1}) via het lineaire transport probleem.
      Bij de 1-adische (gewone) metriek op {0,...,N-1}: dit is de
      Earth Mover Distance (EMD).
  (B) Korter: W_1 op de 3-adische metriek (d_3-adisch(i,j) =
      3^{-k} waar k = diepte van de splitsing). Dit is de
      "natuur lijke" metriek voor de K-L-structuur.
  (C) Is de reeks W_1(mu_k, mu_{k+1}) afnemend? -> gradiëntstroom-
      achtig gedrag.
  (D) Verband met f_2: W_1 per stap zou ~ f_2^k moeten afnemen.

Vereenvoudiging: gebruik de 3-adische metriek want de maat leeft op
residuen mod 3^{k-1}. We coarsen mu_{k+1} naar het groepje mu_k
door te middelen over de 3 sub-blokken.
"""
import numpy as np
from math import log2

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
    return v


def coarsen(v_fine):
    """
    Versimpel v_fine (grootte 3N) naar grootte N door per drietal te sommeren.
    Dit is de push-forward van mu_{k+1} naar het niveau van mu_k.
    """
    N3 = v_fine.size
    assert N3 % 3 == 0
    N  = N3 // 3
    return v_fine[:N] + v_fine[N:2*N] + v_fine[2*N:]


def w1_3adic(mu, nu, k):
    """
    1-Wasserstein afstand met de 3-adische metriek op residuen mod 3^{k-1}.

    d_3adic(i, j) = 3^{-(diepte van de eerste splitsing tussen i en j)}
    In de boom: i en j splitsen op niveau l als i//3^l != j//3^l.

    Vereenvoudiging: we coarsen beide maten hiërarchisch en meten de
    totale variatie op elk niveau (gerelateerd aan W_1 in de ultrametrische zin).

    W_1^{3adic} = Σ_{l=0}^{k-1} 3^{-l} * TV(mu^{(l)}, nu^{(l)})
    waarbij mu^{(l)} = pushforward van mu naar niveau l (grootte 3^l knopen).
    """
    Nl = len(mu)
    assert len(mu) == len(nu) == 3 ** (k - 1), f"{len(mu)} vs {3**(k-1)}"

    tv_sum = 0.0
    m = mu.copy()
    n = nu.copy()
    weight = 1.0

    for l in range(k - 1, -1, -1):
        tv = 0.5 * np.sum(np.abs(m - n))
        tv_sum += weight * tv
        weight *= 3.0
        # Coarsen een niveau
        if len(m) >= 3:
            m = m[:len(m)//3] + m[len(m)//3:2*len(m)//3] + m[2*len(m)//3:]
            n = n[:len(n)//3] + n[len(n)//3:2*len(n)//3] + n[2*len(n)//3:]

    return tv_sum


def kl_divergence(mu, nu):
    """KL-divergentie D_KL(mu || nu) = Σ mu_i log(mu_i/nu_i)."""
    eps = 1e-300
    return float(np.sum(mu * np.log((mu + eps) / (nu + eps))))


print(f"Methode 4: Wasserstein-gradiëntstroom  (lam={LAM})", flush=True)
print("=" * 65, flush=True)

# Bereken consecutieve W1-afstanden en KL-divergenties
print("\n(A) Consecutieve W1 (3-adisch) tussen mu_k en mu_{k+1} geco-arsd",
      flush=True)
print("  k->k+1    W1_3adic(mu_k, mu_{k+1}^coarse)   KL(mu_k || mu_{k+1}^c)",
      flush=True)

prev_w1 = None
for k in range(11, 17):
    v_k    = perron(k)
    v_kp1  = perron(k + 1)
    v_kp1c = coarsen(v_kp1)
    v_kp1c /= v_kp1c.sum()

    # W1 met 3-adische metriek
    w1 = w1_3adic(v_k, v_kp1c, k)
    kl = kl_divergence(v_k, v_kp1c)

    ratio = w1 / prev_w1 if prev_w1 is not None else float('nan')
    print(f"  {k}->{k+1}:    {w1:.6f}   KL={kl:.6f}   "
          f"ratio={ratio:.4f}", flush=True)
    prev_w1 = w1

print("\n(B) TV-afstand per niveau (hiërarchische ontbinding)", flush=True)
print("  k=13 vs k=14 coarsd:", flush=True)
v13  = perron(13)
v14  = perron(14)
v14c = coarsen(v14)
v14c /= v14c.sum()

# Verfijn: TV per 3-adisch niveau
m = v13.copy()
n = v14c.copy()
for l in range(12, -1, -1):
    tv = 0.5 * float(np.sum(np.abs(m - n)))
    print(f"    niveau {l} (grootte {len(m)}): TV = {tv:.6f}", flush=True)
    if len(m) >= 3:
        m = m[:len(m)//3] + m[len(m)//3:2*len(m)//3] + m[2*len(m)//3:]
        n = n[:len(n)//3] + n[len(n)//3:2*len(n)//3] + n[2*len(n)//3:]

print("\n(C) Verband met f2: log(W1) vs k lineair?", flush=True)
w1s = []
for k in range(11, 16):
    v_k    = perron(k)
    v_kp1  = perron(k + 1)
    v_kp1c = coarsen(v_kp1)
    v_kp1c /= v_kp1c.sum()
    w1s.append(w1_3adic(v_k, v_kp1c, k))

log_w1s = np.log(w1s)
ks = np.arange(11, 16, dtype=float)
slope, intercept = np.polyfit(ks, log_w1s, 1)
print(f"  log(W1) = {slope:.4f} * k + {intercept:.4f}", flush=True)
print(f"  Impliciete W1-krimp per stap: exp({slope:.4f}) = {np.exp(slope):.4f}",
      flush=True)
print(f"  Vergelijking f2 ~ 0.90, d_k ~ 0.77", flush=True)

print("\ndone", flush=True)
