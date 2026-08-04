"""
215_eigenvector_predictor.py
============================
Methode 6: Machine learning eigenvector-predictor (sklearn i.p.v. GNN).

Train een Random Forest om v_i^(k) te voorspellen van LOKALE structurele
features van knoop i:
  - Restklasse r_i (0, 1, of 2) -- type van de knoop
  - D3-telling langs het voorouderpad (eerste G stappen omhoog)
  - D1-telling langs het voorouderpad
  - 3-adische valuatie van i (v_3(i) = exponent van 3 in i)
  - "Diepte in de boom" t.o.v. het onderste niveau
  - Aanwezigheid van B3-lift in het voorouderpad

DOEL:
  (A) Hoe goed kan het model v_i voorspellen van lokale features?
      -> R² als maat voor "structurele verklaarbaarheid"
  (B) Feature importance: welke features zijn het meest predictief?
      -> Bevestigt of weerlegt onze analytische inzichten
  (C) Generalisatie: train op k=14, test op k=15
      -> Kunnen we hoger-k eigenvectoren goedkoop benaderen?

Als R² ~ 1: de eigenvector is volledig bepaald door lokale structuur.
Als R² ~ 0: de eigenvector is globaal / holistisch.
Feature importance ~ D3-count: bevestigt Obs 412 (branchingfrequentie).
"""
import numpy as np
from math import log2
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score
from sklearn.inspection import permutation_importance

ALPHA = log2(3.0)
LAM   = 1.70
A     = LAM ** -2.0
B1    = LAM ** (ALPHA - 2.0)
B3    = LAM ** (ALPHA - 1.0)
G     = 6   # voorouder-pad lengte


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
    v /= v.mean()
    return v


def extract_features(k):
    """Extraheer lokale features voor elke knoop i in {0,...,3^{k-1}-1}."""
    N  = 3 ** (k - 1)
    i  = np.arange(N, dtype=np.int64)
    s, r = np.divmod(i, 3)
    Nl = N // 3

    features = {}

    # F1: restklasse (0, 1, 2)
    features['residue'] = r.astype(np.float64)

    # F2: fractie van de index dat deelbaar is door 3 (3-adische structuur)
    def v3(n):
        """3-adische valuatie van n."""
        if n == 0:
            return k
        cnt = 0
        while n % 3 == 0:
            n //= 3
            cnt += 1
        return cnt

    v3_vals = np.array([v3(int(x)) for x in i], dtype=np.float64)
    features['v3_index'] = v3_vals / k  # genormaliseerd

    # F3: D3-telling langs het voorouderpad (G stappen opwaarts in de boom)
    # Voorouder van i op niveau j: i // 3^j
    d3_count = np.zeros(N, dtype=np.float64)
    d1_count = np.zeros(N, dtype=np.float64)
    ancestor = i.copy()
    for step in range(G):
        anc_r = ancestor % 3
        d3_count += (anc_r == 2).astype(np.float64)
        d1_count += (anc_r == 0).astype(np.float64)
        ancestor = ancestor // 3  # omhoog in de boom

    features['d3_anc'] = d3_count / G
    features['d1_anc'] = d1_count / G
    features['d1_minus_d3'] = (d1_count - d3_count) / G

    # F4: niveau van de laagste descendant met r=2 (rijke knopen nabij i)
    # Vereenvoudigd: index mod 9 (de twee niveaus directe structuur)
    features['mod9']  = (i % 9).astype(np.float64) / 9.0
    features['mod27'] = (i % 27).astype(np.float64) / 27.0

    # F5: T4-cyclus-diepte (hoeveel T4-stappen zijn nodig om terug bij i te komen)
    # Vereenvoudigd: i mod 3 (de directe residustructuur)
    features['floor_s'] = (s % Nl).astype(np.float64) / Nl

    # Samenvoegen
    X = np.column_stack([v for v in features.values()])
    feature_names = list(features.keys())
    return X, feature_names


def run_experiment(k_train, k_test=None):
    print(f"\n--- Train op k={k_train}" +
          (f", test op k={k_test}" if k_test else "") + " ---", flush=True)

    X_train, fnames = extract_features(k_train)
    y_train = np.log(perron(k_train))   # log-schaal (beter voor RF)

    # Split voor validatie
    Xtr, Xval, ytr, yval = train_test_split(X_train, y_train,
                                             test_size=0.2, random_state=42)

    rf = RandomForestRegressor(n_estimators=100, max_depth=8,
                                n_jobs=-1, random_state=42)
    rf.fit(Xtr, ytr)

    y_pred_val = rf.predict(Xval)
    r2_val = r2_score(yval, y_pred_val)
    print(f"  R² validatie (k_train): {r2_val:.5f}", flush=True)

    # Feature importance
    imp = rf.feature_importances_
    order = np.argsort(imp)[::-1]
    print(f"  Feature importance:", flush=True)
    for idx in order:
        print(f"    {fnames[idx]:15s}: {imp[idx]:.4f}", flush=True)

    # Generalisatie naar hoger k
    if k_test is not None:
        X_test, _ = extract_features(k_test)
        y_test = np.log(perron(k_test))
        y_pred_test = rf.predict(X_test[:len(y_test)])
        r2_test = r2_score(y_test, y_pred_test)
        print(f"  R² generalisatie (k_test={k_test}): {r2_test:.5f}", flush=True)

    return r2_val


print(f"Methode 6: ML eigenvector-predictor  (lam={LAM})", flush=True)
print("=" * 65, flush=True)

run_experiment(k_train=13, k_test=14)
run_experiment(k_train=14, k_test=15)

print("\ndone", flush=True)
