"""
216_ergodic_3adic.py
====================
Methode 7: Ergodische theorie op de 3-adische solenoïde.

De Collatz-map T: 2N+1 -> (3n+1)/2^{v_2(3n+1)} heeft een uitbreiding
naar de 3-adische gehele getallen Z_3.

Op Z_3: elk element heeft een 3-adische ontwikkeling n = Σ a_i 3^i.
De Collatz-map in residue mod 3^k hangt alleen af van n mod 3^k (bewezen).

Ergodische vraag: als we een willekeurig getal n kiezen en de baan
{n, T(n), T²(n), ...} volgen, convergeert het tijdgemiddelde van
een functie f naar hetzelfde getal (de ruimtegemiddelde E[f]) voor
bijna alle n?

Methode: meet empirisch ergodische gemiddelden.
  - Neem een grote steekproef van startpunten n ~ Uniform[1, M]
  - Volg elke baan T^t(n) voor t = 0..T_max stappen
  - Meet cesaro-gemiddelde: (1/T) Σ_t f(T^t(n))
  - Vergelijk met ruimtegemiddelde E_mu[f] (Perron-maat)
  - Als beide overeenkomen: empirische ergodische aanwijzing

Functies f die we testen:
  (A) f(n) = 1[n ≡ 0 mod 3]   (D1-type indicator)
  (B) f(n) = 1[n ≡ 2 mod 3]   (D3-type indicator)
  (C) f(n) = log(n) / log(M)  (genormaliseerde grootte)
  (D) f(n) = v_2(n+1) / log_2(n)  (halvering-frequentie)

Verwacht: als de Perron-maat de ergodische maat is, dan
E_baan[f] -> E_perron[f] voor bijna alle startpunten.
"""
import numpy as np
from math import log2

ALPHA = log2(3.0)
LAM   = 1.70
A     = LAM ** -2.0
B1    = LAM ** (ALPHA - 2.0)
B3    = LAM ** (ALPHA - 1.0)


def collatz_step(n):
    """Één Collatz-stap (oneven -> (3n+1)/2^k tot oneven)."""
    if n % 2 == 0:
        while n % 2 == 0:
            n //= 2
    else:
        n = 3 * n + 1
        while n % 2 == 0:
            n //= 2
    return n


def collatz_orbit_stats(n0, T_max=200):
    """Volg baan en registreer typefrequenties en groottegemiddelde."""
    n = n0
    d1_cnt = 0
    d3_cnt = 0
    log_sum = 0.0
    halving_sum = 0

    for _ in range(T_max):
        if n == 1:
            break
        r = n % 3
        if r == 0:
            d1_cnt += 1
        elif r == 2:
            d3_cnt += 1

        log_sum += np.log(n)
        # Halverings tellen voor de volgende stap
        m = 3 * n + 1
        hv = 0
        while m % 2 == 0:
            m //= 2
            hv += 1
        halving_sum += hv

        n = collatz_step(n)

    t = d1_cnt + d3_cnt + (T_max - d1_cnt - d3_cnt)
    return d1_cnt, d3_cnt, log_sum, halving_sum, T_max


def perron_type_fractions(k, n_iter=300):
    """Fractie type-0, type-1, type-2 knopen gewogen door Perron-vector."""
    N  = 3 ** (k - 1)
    i  = np.arange(N, dtype=np.int64)
    s, r = np.divmod(i, 3)
    Nl = N // 3
    m0, m1, m2 = (r == 0), (r == 1), (r == 2)

    v = np.ones(N, dtype=np.float64)
    for _ in range(n_iter):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        w  = A * v[T4 := (4 * i + 2) % N]
        w[m2] += B3 * cb[((2 * (i // 3) + 1) % Nl)[m2]]
        w[m0] += B1 * cb[((4 * (i // 3)) % Nl)[m0]]
        v = w / w.max()
    v /= v.sum()

    f0 = float(v[m0].sum())
    f1 = float(v[m1].sum())
    f2 = float(v[m2].sum())
    return f0, f1, f2


print(f"Methode 7: Ergodische theorie 3-adisch  (lam={LAM})", flush=True)
print("=" * 65, flush=True)

# Perron-maat typefrequenties (k=13)
k_ref = 13
N_ref = 3 ** (k_ref - 1)
i_ref = np.arange(N_ref, dtype=np.int64)
T4_ref = (4 * i_ref + 2) % N_ref
s_ref, r_ref = np.divmod(i_ref, 3)
Nl_ref = N_ref // 3
m0_ref, m1_ref, m2_ref = (r_ref == 0), (r_ref == 1), (r_ref == 2)
v_ref = np.ones(N_ref, dtype=np.float64)
for _ in range(300):
    cb_r = np.minimum(np.minimum(v_ref[:Nl_ref], v_ref[Nl_ref:2*Nl_ref]),
                       v_ref[2*Nl_ref:])
    w_r  = A * v_ref[T4_ref]
    w_r[m2_ref] += B3 * cb_r[((2 * s_ref + 1) % Nl_ref)[m2_ref]]
    w_r[m0_ref] += B1 * cb_r[((4 * s_ref) % Nl_ref)[m0_ref]]
    v_ref = w_r / w_r.max()
v_ref /= v_ref.sum()

p_d1 = float(v_ref[m0_ref].sum())  # type 0 = D1 in omgekeerde richting
p_d3 = float(v_ref[m2_ref].sum())  # type 2 = D3

print(f"\nPerron-maat typefrequenties (k={k_ref}):", flush=True)
print(f"  P(type-0/D1) = {p_d1:.5f}", flush=True)
print(f"  P(type-1)    = {float(v_ref[m1_ref].sum()):.5f}", flush=True)
print(f"  P(type-2/D3) = {p_d3:.5f}", flush=True)

# Empirische baangemiddelden
print(f"\nEmpirische baangemiddelden (1000 willekeurige startpunten, T=200):",
      flush=True)

rng = np.random.default_rng(42)
M   = 10**6
starts = rng.integers(1, M, size=1000)
starts = starts | 1  # maak oneven

d1_fracs = []
d3_fracs = []
halving_fracs = []

for n0 in starts:
    d1, d3, logs, hv, T = collatz_orbit_stats(int(n0), T_max=200)
    total = d1 + d3 + (T - d1 - d3)
    if total > 10:
        d1_fracs.append(d1 / T)
        d3_fracs.append(d3 / T)
        halving_fracs.append(hv / T)

print(f"  E_baan[type-0/D1] = {np.mean(d1_fracs):.5f}  "
      f"std={np.std(d1_fracs):.5f}  Perron={p_d1:.5f}", flush=True)
print(f"  E_baan[type-2/D3] = {np.mean(d3_fracs):.5f}  "
      f"std={np.std(d3_fracs):.5f}  Perron={p_d3:.5f}", flush=True)
print(f"  E_baan[halvering/stap] = {np.mean(halving_fracs):.5f}  "
      f"std={np.std(halving_fracs):.5f}", flush=True)
print(f"  Theorie halvering/stap ~ log2(3) / 2 = {log2(3)/2:.5f}", flush=True)

# Test ergodische convergentie: variantie van baangemiddelden
print(f"\nErgodische variantie (variantie van baangemiddelden):", flush=True)
print(f"  Var(D3_frac per baan) = {np.var(d3_fracs):.6f}", flush=True)
print(f"  Als ergodisch: var -> 0 voor T->inf", flush=True)

# Vergelijk korte vs lange banen
d3_short = []
d3_long  = []
for n0 in starts[:200]:
    d1s, d3s, _, _, T = collatz_orbit_stats(int(n0), T_max=50)
    d1l, d3l, _, _, T = collatz_orbit_stats(int(n0), T_max=200)
    d3_short.append(d3s / 50)
    d3_long.append(d3l / 200)

print(f"  Var(D3, T=50)  = {np.var(d3_short):.6f}", flush=True)
print(f"  Var(D3, T=200) = {np.var(d3_long):.6f}", flush=True)
print(f"  Ratio: {np.var(d3_long)/np.var(d3_short):.4f} "
      f"(ergodisch -> 50/200 = 0.25)", flush=True)

print("\ndone", flush=True)
