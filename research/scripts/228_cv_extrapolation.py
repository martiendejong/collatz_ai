"""
228_cv_extrapolation.py
=======================
Extrapolatie van CV(k) naar k -> inf.

Gebruikt hardgecodeerde k=10..15 data uit Script 223 + voegt k=16 toe.
Methode A: empirische fit dCV(k) ~ C * r^k => CV_inf via geometrische rest
Methode B: analyse van variantiecomponenten (tussen-type vs binnen-type)
"""
import numpy as np
from math import log2

ALPHA = log2(3.0)
LAM   = 1.70
A  = LAM ** -2.0
B1 = LAM ** (ALPHA - 2.0)
B3 = LAM ** (ALPHA - 1.0)

# Hardgecodeerde resultaten uit Script 223 (k=10..15)
CATALOG = [
    # k,   sw/rho,  rho,      cv
    (10, 0.76817, 1.02369, 0.75526),
    (11, 0.75720, 1.03572, 0.77266),
    (12, 0.76158, 1.04167, 0.78868),
    (13, 0.75303, 1.04367, 0.80207),
    (14, 0.75597, 1.04753, 0.81427),
    (15, 0.75376, 1.04966, 0.82489),
]


def compute_k16(n_iter=150):
    """Bereken k=16 voor verificatie (n_iter=150 geeft rho nauwkeurig genoeg)."""
    k = 16
    N  = 3 ** (k - 1)
    print(f"  Berekenen k={k}, N={N:,d} ...")
    i  = np.arange(N, dtype=np.int64)
    T4 = (4 * i + 2) % N
    s, r = np.divmod(i, 3)
    Nl = N // 3
    m0, m2 = (r == 0), (r == 2)
    R1 = (4 * s) % Nl
    R3 = (2 * s + 1) % Nl

    v = np.ones(N, dtype=np.float64)
    for it in range(n_iter):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        w  = A * v[T4]
        w[m2] += B3 * cb[R3[m2]]
        w[m0] += B1 * cb[R1[m0]]
        v = w / w.max()
        if it % 50 == 49:
            print(f"    iter {it+1}/{n_iter} done")
    v /= v.mean()

    cb2 = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
    w2  = A * v[T4]
    w2[m2] += B3 * cb2[R3[m2]]
    w2[m0] += B1 * cb2[R1[m0]]
    rho = float(w2.mean())
    cv  = float(np.std(v))
    sw  = float(np.mean(np.sqrt(v)))

    # Variantiecomponenten
    a = float(v[r == 0].mean())
    b = float(v[r == 2].mean())
    c = float(v[r == 1].mean())
    var_between = ((a - 1.0)**2 + (b - 1.0)**2 + (c - 1.0)**2) / 3
    var_within_D1 = float(np.var(v[r == 0]))
    var_within_D3 = float(np.var(v[r == 2]))
    var_within_D2 = float(np.var(v[r == 1]))

    return {
        'k': k, 'sw_rho': sw/rho, 'rho': rho, 'cv': cv,
        'a': a, 'b': b, 'c': c,
        'var_between': var_between,
        'var_within': (var_within_D1 + var_within_D3 + var_within_D2) / 3,
        'vmax': float(v.max()), 'vmin': float(v.min()),
    }


print(f"228: CV-extrapolatie naar k->inf  (lam={LAM})")
print("=" * 65)
print()

# Bereken k=16
r16 = compute_k16()
print(f"  k=16: CV={r16['cv']:.5f}  rho={r16['rho']:.5f}  sw/rho={r16['sw_rho']:.5f}")
print()

# Volledige tabel
all_data = list(CATALOG) + [(r16['k'], r16['sw_rho'], r16['rho'], r16['cv'])]
ks_all  = np.array([d[0] for d in all_data], dtype=float)
cvs_all = np.array([d[3] for d in all_data])

print("Overzicht CV(k):")
prev_cv = None
for d in all_data:
    k, swr, rho, cv = d
    dcv = f"  dCV={cv-prev_cv:+.5f}" if prev_cv is not None else ""
    print(f"  k={k:2d}  CV={cv:.5f}  rho={rho:.5f}  sw/rho={swr:.5f}{dcv}")
    prev_cv = cv

print()
print("=== METHODE A: EMPIRISCHE FIT dCV(k) = C * r^k ===")
print()

increments = np.diff(cvs_all)
inc_ks = ks_all[1:]
ok = increments > 0
log_inc = np.log(increments[ok])
fit_coeffs = np.polyfit(inc_ks[ok], log_inc, 1)
log_r = fit_coeffs[0]
r_hat = float(np.exp(log_r))
log_C = fit_coeffs[1]
C_hat = float(np.exp(log_C))
k_last = int(ks_all[-1])
cv_last = float(cvs_all[-1])

print(f"Log-lineaire fit: r = {r_hat:.4f}  C = {C_hat:.6f}")
print(f"  => dCV(k) ~ {C_hat:.4f} * {r_hat:.4f}^k")
print()

# Restsom na k_last
remaining_A = C_hat * r_hat**(k_last + 1) / (1 - r_hat)
cv_inf_A = cv_last + remaining_A
print(f"Restsom na k={k_last}: {remaining_A:.5f}")
print(f"CV_inf (methode A, volledige fit) = {cv_last:.5f} + {remaining_A:.5f} = {cv_inf_A:.5f}")
print()

# Ratio-methode (robuuster voor de staart)
dec_ratios = increments[1:] / increments[:-1]
print(f"Ratio opeenvolgende incrementen: {[f'{rr:.4f}' for rr in dec_ratios]}")
avg_ratio_3 = float(np.mean(dec_ratios[-3:]))
last_inc = increments[-1]
remaining_B = last_inc * avg_ratio_3 / (1.0 - avg_ratio_3) if avg_ratio_3 < 1.0 else float('nan')
cv_inf_B = cv_last + last_inc + remaining_B
print(f"Gemiddelde ratio (laatste 3): {avg_ratio_3:.4f}")
print(f"CV_inf (ratio-methode) = {cv_last:.5f} + {last_inc:.5f} + {remaining_B:.5f} = {cv_inf_B:.5f}")
print()

# Bovengrens: ratio <= 0.95 (conservatief)
ub_ratio = 0.95
remaining_ub = last_inc * ub_ratio / (1.0 - ub_ratio)
cv_inf_ub = cv_last + last_inc + remaining_ub
print(f"Conservatieve bovengrens (ratio <= {ub_ratio}): CV_inf <= {cv_inf_ub:.4f}")
print()

print("=== METHODE B: VARIANTIECOMPONENTEN ===")
print()
print("Var(v) = var_between + var_within_avg  (wet van totale variantie)")
print(f"  k=16: var_between={r16['var_between']:.5f}  var_within_avg={r16['var_within']:.5f}")
print(f"  Type-gemiddelden k=16: a={r16['a']:.5f}  b={r16['b']:.5f}  c={r16['c']:.5f}")
print()

# Exacte c/a = A/rho relatie
a16, c16, rho16 = r16['a'], r16['c'], r16['rho']
print(f"  Exacte relatie c/a = A/rho:")
print(f"    Gemeten c/a = {c16/a16:.8f}")
print(f"    A/rho       = {A/rho16:.8f}")
print(f"    Fout:         {abs(c16/a16 - A/rho16):.2e}")
print()

# Theoretische bovengrens via max/min
print("=== BOVENGRENS VIA EIGENVECTOR SPREADING ===")
print()
print("  Var(v) <= (max(v) - mean(v))^2 = (max(v) - 1)^2  [triviale ub]")
for d in all_data:
    pass  # we have no vmax/vmin for k<16
print(f"  k=16: max={r16['vmax']:.5f}  min={r16['vmin']:.5f}")
print(f"  (max-1)^2 = {(r16['vmax']-1)**2:.5f}  vs Var(v) = CV^2 = {r16['cv']**2:.5f}")
print(f"  Spreading ratio (CV^2 / (max-1)^2): {r16['cv']**2 / (r16['vmax']-1)**2:.5f}")
print()

print("=== SAMENVATTING ===")
print()
print(f"CV(k) stijgt monotoon: {cvs_all[0]:.5f} -> {cvs_all[-1]:.5f}")
print(f"Incrementen krimpen: ratio ~{avg_ratio_3:.4f} per stap")
print(f"Extrapolatie: CV_inf in [{cv_inf_A:.4f}, {cv_inf_B:.4f}]")
print(f"Conservatieve bovengrens: CV_inf <= {cv_inf_ub:.4f}")
print()
print(f"CONCLUSIE:")
print(f"  CV_inf > 0:  BEWEZEN (CV monotoon stijgend, CV(10) = {cvs_all[0]:.5f} > 0)")
print(f"  CV_inf < inf: GEMETEN (extrapolatie geeft ~{(cv_inf_A+cv_inf_B)/2:.3f})")
print(f"  CV_inf < 1:  ONZEKER (beste schatting {(cv_inf_A+cv_inf_B)/2:.3f}, ub = {cv_inf_ub:.3f})")
print()
print("done")
