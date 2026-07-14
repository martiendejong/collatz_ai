"""E24: cross-analysis of all program data — consistency web, records law, unification checks."""
import sys, os, json, math
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "lib"))
sys.stdout.reconfigure(encoding="utf-8")
from collatz_lib import *
out = {}

# ===== A. the consistency web: every constant from ONE source (geometric(1/2) + log2 3) =====
L = math.log2(3)
web = {
 "syracuse_drift": dict(derived=round(L - 2, 4), measured=-0.415),
 "macro_drift(=E[k]x)": dict(derived=round(2 * (L - 2), 4), measured=-0.83),
 "ternlen_drift(xlog3_2)": dict(derived=round((L - 2) / L, 4), measured=-0.26),
 "ladder_rung(P(k=2)P(w=1))": dict(derived=0.125, measured=0.125),
 "excursion_tail(martingale)": dict(derived=-1.0, measured=-0.978),
 "alpha_star": dict(derived=round((1 + L) / 2, 5), measured=1.3),
 "spike_slope": dict(derived=round((L - 1) / 2, 5), measured=0.29),
 "needle_endpoint_KL": dict(derived=round((1/L)*math.log2(2/L)+(1-1/L)*math.log2(2*(1-1/L)), 4), measured=0.05),
 "streak_cont(~P(nondec))": dict(derived=0.70, measured=0.68),
}
out["consistency_web"] = web

# ===== B. unification: 3n-1 on positives == 3n+1 on negatives (conjugacy n <-> -n) =====
ok = True
import random
random.seed(61)
for _ in range(10000):
    n = random.randrange(1, 1 << 40)
    a = 3 * n - 1 if n % 2 else n // 2          # 3n-1 map on positive n
    b = -(3 * (-n) + 1) if n % 2 else -((-n) // 2)  # mirrored Collatz on -n
    if a != b: ok = False; break
out["conjugacy_3nm1_eq_neg_collatz"] = ok
out["unification"] = "3n-1 cycles {1,5,17} == census realizations {-1,-5,-17}: the control system IS the negative half; Sign Theorem explains why the climbers landed there"

# ===== C. records law: t_rec vs log2(n) linear fit + prediction =====
recs = [(27,59),(703,81),(10087,105),(35655,135),(270271,164),(362343,165),
        (381727,173),(626331,176),(1027431,183),(1126015,224)]
xs = [math.log2(n) for n,_ in recs]; ys = [t for _,t in recs]
mx, my = sum(xs)/len(xs), sum(ys)/len(ys)
slope = sum((x-mx)*(y-my) for x,y in zip(xs,ys)) / sum((x-mx)**2 for x in xs)
icept = my - slope*mx
out["records_law"] = dict(slope_per_bit=round(slope,2), intercept=round(icept,1),
    theory_slope="1/0.080 = 12.5 (tail exponent inverse)",
    prediction_t300=f"first stopping-time >= 300 near n ~ 2^{round((300-icept)/slope,1)} ~ {2**((300-icept)/slope):.2e}")

# ===== D. champions' coordinates: where is their capital? =====
champs = {}
for n,t in recs:
    a,k = coords(n)
    champs[n] = dict(k=k, a_bits=a.bit_length(), stop=t)
out["champion_coords"] = champs
out["champion_reading"] = "k in 2..7 (elevated vs mean 2 but modest): record capital is distributed along the orbit, not stored upfront -- consistent with F6/F11"

# ===== E. one new pattern probe: do the three integer cycles' a-values (1,-1,-5) relate to unit/Catalan? =====
out["cycle_a_values"] = {"n=1": 1, "n=-5": -1, "n=-17": -1, "n=-41(phase)": -5}
out["reading"] = "all realized a-values are the divisors of the Catalan-gap era: +-1 (free) and -5; no large-a cycle exists in range -- integrality forces tiny families"
print(json.dumps(out, indent=1))
