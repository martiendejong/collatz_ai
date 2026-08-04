"""
226_task4_maintenance.py
========================
Verificatie van Taak 4: de getilte onderhoudsfactor is < env/(2/3) = 1.10.

Taak 4 stelt: voor dominantie-ketens van lengte g geldt dat de flow-massa
per stap (de "getilte onderhoudsfactor") strikt kleiner is dan env/(2/3).

Ingredienten (allemaal al beschikbaar):
  (A) Freshness-lemma: telmassa per ketensstap <= (2/3) [BEWEZEN]
  (B) C_tilt < 1 voor het volledige veld [GEMETEN in script 225, Obs 430]
  (C) Uit A+B: flowmassa per stap <= (2/3) * C_tilt < (2/3) * 1 = 0.667

Onderhoudsfactor = r(g) / (2/3), gemeten 0.60-0.78 (scripts 195, 197).
env/(2/3) = 0.731..0.735 / (2/3) = 1.097..1.103.

Doel: maat r(g) direct (als W{keten>=g+1}/W{keten>=g}) en verifieer
dat r(g) < env voor alle (k, eps, g).

Definitie:
  G(m) = F(4m+2) - F(m) = log2(v(4m+2)) - log2(v(m))
  Dominantie op niveau g: G(m) <= -t_0(eps) voor g opeenvolgende stappen
  (met t_0(eps) = -log2(eps * lam^2 * rho))

  W{keten >= g} = E_W[1{G_1 <= -t0, ..., G_g <= -t0}]
               = (1/N) * sum_{keten-startpunten} v(m) * 1{keten>=g}
               (want W = flow-maat = v/E[v] en E[v]=1)

Meting: voor k=11,13,15 en eps=0.05,0.10:
  - Bereken chain_mass_g = (1/N) * sum v(m) * 1{keten>=g}
  - r(g) = chain_mass_{g+1} / chain_mass_g
  - Vergelijk r(g) met env en env/(2/3)
"""
import numpy as np
from math import log2, log

ALPHA = log2(3.0)
LAM   = 1.70
A  = LAM ** -2.0
B1 = LAM ** (ALPHA - 2.0)
B3 = LAM ** (ALPHA - 1.0)
ENV_ANALYTIC = (B1 + B3) / 3  # rho_1 = 0.7221; env = (B1+B3)/(3*rho) ~ 0.733


def run_analysis(k, eps_list=(0.05, 0.10), n_iter=400, g_max=8):
    N  = 3 ** (k - 1)
    i  = np.arange(N, dtype=np.int64)
    T4 = (4 * i + 2) % N
    s, r = np.divmod(i, 3)
    Nl = N // 3
    m0, m2 = (r == 0), (r == 2)
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
    env = (B1 + B3) / (3 * rho)

    # F = log2(v)
    F = np.log2(v)

    print(f"\n  k={k}, N={N}, rho={rho:.5f}, env={env:.5f}, env/(2/3)={env/(2/3):.5f}")

    for eps in eps_list:
        t0 = -log2(eps * LAM**2 * rho)   # drempelwaarde

        # G(m) = F(4m+2) - F(m) voor m in het eerste derde (0..Nl-1)
        # m in Nl: backbone m -> 4m+2 mod N (maar dit is T4 van het hele geval)
        # Voor s in {0..Nl-1}: G(s) = F(T4_even(s)) - F(s)
        # T4_even: s -> (4s+2) mod Nl? Nee: T4(i) voor i=3s+r geeft (4i+2) mod N.
        # Backbone: m zit in een r=0 of r=2 positie; de feed map gaat naar Nl.
        # Voor het dominantiebegrip: G = F(backbone_succes) - F(current_pos)
        # Eenvoudigst: gebruik T4 zelf als de "step" in de keten.
        # Keten: m -> T4(m) -> T4(T4(m)) -> ...
        # G_step(m) = F(T4(m)) - F(m)
        G = F[T4] - F  # G(m) = log2(v[T4(m)]) - log2(v[m])

        # Dominantie: G(m) <= -t0
        dom = (G <= -t0)

        # Ketenmassa: W{keten >= g} = (1/N) * sum_{m: dom_1,...,dom_g} v(m)
        chain_mask = np.ones(N, dtype=bool)
        chain_masses = []
        for g in range(g_max + 1):
            mass = float(np.sum(v[chain_mask]) / N)
            chain_masses.append(mass)
            if mass < 1e-12: break
            # Update keten: next stap vereist ook G(T4^g(m)) <= -t0
            # Maar we kijken naar ketens die STARTEN bij m en g stappen doen.
            # Eenvoudiger: iteratief bijhouden welke start-punten keten>=g hebben.
            chain_mask = chain_mask & dom
            # Schuif: kijk of T4(m) in chain_mask zit voor de volgende stap
            next_mask = np.zeros(N, dtype=bool)
            next_mask[T4[chain_mask]] = True   # markeringen van de volgende posities
            # Nee: dit klopt niet. chain_mask geeft startpunten waarvan de keten >= g stappen gaat.
            # Herdefinieer: chain_mask[m] = True als G(m) <= -t0 EN keten al g stappen doorlopen.
            # Begin: chain_mask = alle m (keten >= 0).
            # Na stap 1: chain_mask[m] = dom[m] (G(m) <= -t0).
            # Na stap 2: chain_mask[m] = dom[m] AND dom[T4(m)].
            # Etc.
            chain_mask = chain_mask & dom[T4]  # Wacht, dit werkt niet goed.
            # Correcte logica: bewaar "heeft_keten_g" als boolean array over ALLE m.
            # heeft_keten_0 = True voor alle m (triviale keten)
            # heeft_keten_{g+1}[m] = heeft_keten_g[m] AND dom[m]
            # Hmm maar dan is chain_mass_g = E_W[heeft_keten_g] maar met dom op ANDERE plaatsen.
            break  # de huidige aanpak is fout, herstart

        # Correcte ketenmassa-berekening:
        chain_masses = []
        has_chain_g = np.ones(N, dtype=bool)  # keten >= 0 (altijd)
        for g in range(g_max + 1):
            mass = float(np.mean(v[has_chain_g]) * has_chain_g.mean())
            # = sum_m v[m] * has_chain_g[m] / N
            mass2 = float(np.sum(v * has_chain_g) / N)
            chain_masses.append(mass2)
            if mass2 < 1e-15: break
            # Keten >= g+1: heeft_keten_g EN dom op de g+1-de positie
            # Positie g+1 voor startpunt m is T4^{g+1}(m).
            # Iteratief: pos_g = T4^g(m). Begin pos_0 = m.
            # Eenvoudigste aanpak: bijhouden welke startpunten nog >= g stappen kunnen gaan.
            # Na stap: kijk of de huidige positie (na g stappen) dom heeft.
            # But chain_mask is over startpunten... we need the CURRENT position after g steps.
            # Let's just track the "current position" for each still-surviving chain.

        # Helemaal opnieuw: bijhoud voor elk startpunt de positie na g stappen.
        chain_masses = []
        current_pos = np.arange(N, dtype=np.int64)  # beginpositie = zichzelf
        surviving = np.ones(N, dtype=bool)  # alle ketens overleven bij g=0
        for g in range(g_max + 1):
            # Massa bij g: som van v[startpunt] voor alle overlevende ketens
            mass_g = float(np.sum(v[surviving]) / N)
            chain_masses.append(mass_g)
            if mass_g < 1e-15: break
            # Update: overleving vereist G(current_pos) <= -t0
            # current_pos = T4^g(startpunt)
            surviving = surviving & dom[current_pos]
            current_pos = T4[current_pos]

        # Bereken r(g) = chain_mass_{g+1} / chain_mass_g
        ratios = []
        for g in range(len(chain_masses) - 1):
            if chain_masses[g] > 1e-15:
                r_g = chain_masses[g+1] / chain_masses[g]
                ratios.append(r_g)
            else:
                break

        # Verificatie: r(g) vs env en env/(2/3)
        max_r = max(ratios) if ratios else float('nan')
        print(f"    eps={eps:.2f}  t0={t0:.3f}  chain_masses={[f'{m:.4f}' for m in chain_masses[:6]]}")
        print(f"    r(g)={[f'{r:.4f}' for r in ratios[:6]]}  max_r={max_r:.4f}  "
              f"vs env={env:.4f}  env/(2/3)={env/(2/3):.4f}  margin={env/max_r:.3f}x")

    return env


print(f"226: Taak 4 verificatie getilte onderhoudsfactor  (lam={LAM})")
print("=" * 70)
print()
print(f"Analytisch: env = (B1+B3)/3 / rho = {ENV_ANALYTIC:.5f}/rho")
print(f"Bewijs-route: r(g) < env omdat C_tilt < 1 (T3) en telratio <= 2/3 (Freshness)")
print(f"Grens: tilted_maintenance = r(g)/(2/3) < env/(2/3)")
print()

for k in (11, 13, 15):
    run_analysis(k)

print()
print("=== TAAK 4 STATUS ===")
print("""
Als de metingen r(g) < env tonen voor alle (k,eps,g), dan is Taak 4 GESLOTEN:
  flow_maintenance = r(g) <= C_tilt * (2/3) <= 1 * (2/3) < env (want (2/3) < env ~ 0.733)

De T3-resultaat (C_tilt < 1) geeft hier DIRECT de gewenste bound:
  r(g) <= C_tilt * count_ratio <= 1 * env^* < env

waarbij count_ratio <= env^* (het telling-per-stap gemiddelde) via Freshness.
""")
print("done")
