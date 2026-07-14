"""E6b fix: high-precision continued fraction of log2(3) and the cycle-length bound."""
import sys, json
sys.stdout.reconfigure(encoding="utf-8")
from decimal import Decimal, getcontext
from fractions import Fraction

getcontext().prec = 90
lg = Decimal(3).ln() / Decimal(2).ln()   # log2(3), 90 digits
out = {}

# continued fraction expansion
x = lg
cf = []
for _ in range(40):
    ai = int(x)
    cf.append(ai)
    fr = x - ai
    if fr == 0: break
    x = 1 / fr
out["cf_terms"] = cf[:20]

# convergents p/q
ps, qs = [0, 1], [1, 0]
convs = []
lgF = Fraction(str(lg))
for ai in cf:
    ps.append(ai * ps[-1] + ps[-2])
    qs.append(ai * qs[-1] + qs[-2])
    p, q = ps[-1], qs[-1]
    err = abs(lgF - Fraction(p, q))
    convs.append((p, q, float(err)))
out["convergents"] = [(p, q, f"{e:.3e}") for p, q, e in convs[:18]]

# cycle constraint: 2^D = 3^K * prod(1+1/(3 n_i)), n_i > 2^71 (verification floor), K odd steps.
# => 0 < D - K*log2(3) < K * log2(e) / (3*2^71)  =>  |log2(3) - D/K| < 4.81e-22  (per unit)
eps = float(Fraction(1) * 1)
import math
eps = math.log2(math.e) / (3 * 2 ** 71)
out["required_approx_error"] = f"{eps:.3e}"
# legendre: only convergents can approximate better than 1/(2 q^2); error of convergent q_i ~ 1/(q_i q_{i+1})
# need 1/(q_i * q_{i+1}) < eps  -> find first convergent where the NEXT denominator makes it possible
feasible = None
for i in range(len(convs) - 1):
    p, q, e = convs[i]
    if e < eps:
        feasible = (p, q)
        break
out["first_convergent_with_small_enough_error"] = feasible
# K must be >= that q (or a combination); min cycle length = K + D ~ K(1+log2 3)
if feasible:
    q = feasible[1]
    out["min_odd_steps_K"] = q
    out["min_cycle_length"] = int(q * (1 + float(lg)))
print(json.dumps(out, indent=1))
