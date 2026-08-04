"""
211_persistent_homology.py
==========================
Methode 2: Vereenvoudigde aanhoudende homologie op het Collatz-eigenruimtelandschap.

Zonder gudhi: implementeer zelf de Vietoris-Rips filtratie op kleine schaal.
H0 = samengevoegde componenten (drempelwaarde van eigenvector-hoogte)
H1 = lussen in de graaf (potentiële "bijna-cycli")

Aanpak:
- Bouw de Collatz-boom als gerichte graaf op niveau k (N=3^{k-1} knopen)
- Gebruik v^(k) als hoogtef unctie op de knopen
- Bouw subniveaufiltratie: voeg knopen toe in volgorde van stijgende v_i
- Bijhoud: wanneer verschijnen nieuwe samengevoegde componenten (H0-bars)?
- Wanneer sluit een lus (H1-generator)?
- Een H1-generator die lang "leeft" (grote dood-geboorte-ratio) = topologische bijna-cyclus

Vereenvoudiging: gebruik alleen de Collatz-boom-richtingen als edges
(niet de volledige Rips-complex), want de boom-topologie is het relevante object.
"""
import numpy as np
from math import log2

ALPHA = log2(3.0)
LAM   = 1.70
A     = LAM ** -2.0
B1    = LAM ** (ALPHA - 2.0)
B3    = LAM ** (ALPHA - 1.0)


def build_and_perron(k, n_iter=300):
    N  = 3 ** (k - 1)
    i  = np.arange(N, dtype=np.int64)
    T4 = (4 * i + 2) % N
    s, r = np.divmod(i, 3)
    Nl = N // 3
    R1 = (4 * s) % Nl
    R3 = (2 * s + 1) % Nl
    m0, m2 = (r == 0), (r == 2)

    v = np.ones(N, dtype=np.float64)
    for _ in range(n_iter):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        w  = A * v[T4]
        w[m2] += B3 * cb[R3[m2]]
        w[m0] += B1 * cb[R1[m0]]
        v = w / w.max()
    v /= v.mean()
    return v, T4, R1, R3, m0, m2, Nl, N


def union_find(n):
    parent = np.arange(n, dtype=np.int64)
    rank   = np.zeros(n, dtype=np.int64)

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra == rb:
            return False, -1, -1   # lus gevonden
        if rank[ra] < rank[rb]:
            ra, rb = rb, ra
        parent[rb] = ra
        if rank[ra] == rank[rb]:
            rank[ra] += 1
        return True, ra, rb         # samengevoegd

    return find, union


def persistent_homology(k):
    """
    H0 en H1 van de Collatz subniveaufiltratie.
    Edges = de K-L-operatorverbindingen (T4, lift-targets).
    Knopen worden toegevoegd op volgorde van v_i.
    """
    v, T4, R1, R3, m0, m2, Nl, N = build_and_perron(k)

    # Bouw randlijst: alle operatorverbindingen
    edges = []
    i_all = np.arange(N, dtype=np.int64)
    # Walk-edge: i -> T4(i)
    for i in range(N):
        edges.append((i, int(T4[i])))
    # Lift-edges voor type 2
    idx2 = np.where(m2)[0]
    s2   = idx2 // 3
    for idx in idx2:
        s = idx // 3
        for j in range(3):
            edges.append((idx, int(R3[idx] + j * Nl)))
    # Lift-edges voor type 0
    idx0 = np.where(m0)[0]
    for idx in idx0:
        s = idx // 3
        for j in range(3):
            edges.append((idx, int(R1[idx] + j * Nl)))

    edges = list(set(edges))  # dedup

    # Sorteervolgorde: edge birth = max van v-waarden van de twee eindpunten
    # Knoop birth = v_i (toevoegmoment)
    node_birth = v.copy()

    # Sorteer knopen op geboortetijd
    node_order = np.argsort(node_birth)
    node_time  = np.empty(N)
    node_time[node_order] = np.arange(N)  # rang als tijdstempel

    # Edge birth = max(node_birth[u], node_birth[v])
    edge_births = [(max(node_birth[u], node_birth[v]), u, v) for u, v in edges]
    edge_births.sort()

    # Filtratie: voeg knopen en edges toe in volgorde
    find, union = union_find(N)
    active = np.zeros(N, dtype=bool)

    # Interleave knopen en edges gesorteerd op drempelwaarde
    events = []
    for i in range(N):
        events.append((node_birth[i], 0, i, -1))  # knoop
    for b, u, v_node in edge_births:
        events.append((b, 1, u, v_node))  # edge
    events.sort()

    h0_bars = []   # (geboorte, dood) voor H0
    h1_bars = []   # (geboorte, dood) voor H1

    component_birth = {}  # component root -> geboortetijd

    for val, kind, u, v_node in events:
        if kind == 0:
            # Knoop i=u wordt toegevoegd
            active[u] = True
            r = find(u)
            component_birth[r] = val
        else:
            # Edge (u, v_node) wordt toegevoegd
            if not (active[u] and active[v_node]):
                continue
            merged, ra, rb = union(u, v_node)
            if merged:
                # H0: twee componenten worden één; kleinere sterft
                # Geboorte van stervende component
                r_old = rb  # rb is samengevoegd in ra
                b_old = component_birth.get(r_old, val)
                b_new = component_birth.get(ra, val)
                # Stervende component had geboorte min(b_old, b_new)
                dying_birth = min(b_old, b_new)
                h0_bars.append((dying_birth, val))
                # Update birth van overlevende
                component_birth[ra] = min(b_old, b_new)
                if r_old in component_birth:
                    del component_birth[r_old]
            else:
                # Lus gesloten: H1-generator geboren
                h1_bars.append((val, np.inf))  # dood = inf tenzij we boundary toevoegen

    # Overlevende H0-componenten
    n_components = len(component_birth)

    # Analyse H1: tel generators en hun geboortetijden
    # (in onze benadering sterven H1-generators nooit = de lussen zijn persistent)
    h1_births = [b for b, d in h1_bars]

    return h0_bars, h1_bars, n_components, v, node_birth


print(f"Methode 2: Aanhoudende homologie  (lam={LAM})", flush=True)
print("=" * 65, flush=True)

for k in (6, 7, 8):
    print(f"\n--- k={k}  (N={3**(k-1):,}) ---", flush=True)
    h0_bars, h1_bars, n_comp, v, nb = persistent_homology(k)

    # H0 statistieken
    h0_lifetimes = [d - b for b, d in h0_bars if d < np.inf]
    print(f"  H0 bars (samengevoegde comp.): {len(h0_bars)}", flush=True)
    print(f"  H0 gem. levensduur: {np.mean(h0_lifetimes):.5f}  "
          f"max: {np.max(h0_lifetimes):.5f}", flush=True)
    print(f"  Overlevende componenten (k->inf): {n_comp}", flush=True)

    # H1 statistieken
    print(f"  H1 generators (lussen): {len(h1_bars)}", flush=True)

    if h1_bars:
        h1_birth_vals = [b for b, d in h1_bars]
        print(f"  H1 geboorte: min={min(h1_birth_vals):.5f}  "
              f"max={max(h1_birth_vals):.5f}  "
              f"gem={np.mean(h1_birth_vals):.5f}", flush=True)

        # Welk fractie van de H1-generatoren wordt geboren in het onderste kwartiel?
        q25 = np.quantile(nb, 0.25)
        frac_low = np.mean(np.array(h1_birth_vals) < q25)
        print(f"  H1 geboren in onderste kwartiel (bijna-cycli bij lage v): "
              f"{frac_low:.4f}", flush=True)

    # Controleer: H0 = N-1 (boom) + extra vanwege de extra edges?
    print(f"  Verwacht H0 voor boom: {3**(k-1) - 1}", flush=True)

print("\ndone", flush=True)
