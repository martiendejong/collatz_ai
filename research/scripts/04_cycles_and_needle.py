"""E6: cycle equations (circuit eq, convergents of log2 3, length quantization).
E7: never-drop needle DP + survivor structure."""
import sys, os, json, math
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "lib"))
sys.stdout.reconfigure(encoding="utf-8")
from collatz_lib import *
from fractions import Fraction
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

FIG = os.path.join(os.path.dirname(__file__), "..", "figures")
out = {}

# ============ E6a: 1-circuit equation a(2^(k+w)-3^k) = 2^w-1, wide search ============
sols = []
for k in range(1, 400):
    p3 = 3 ** k
    for w in range(1, 400):
        d = (1 << (k + w)) - p3
        if d <= 0: continue
        if ((1 << w) - 1) % d == 0:
            a = ((1 << w) - 1) // d
            if a > 0 and a % 2 == 1:
                sols.append((a, k, w))
out["circuit_solutions_k_w_lt_400"] = sols  # expect only (1,1,1) - Steiner's theorem

# ============ E6b: continued fraction convergents of log2(3): allowed cycle shapes ============
# any cycle: 2^D / 3^K = prod(1 + 1/(3 n_i)) with n_i > 2^71  =>
# 0 < D log2 - K log3 < sum 1/(3 n_i) * (1/ln2 adj) < K / (3 * 2^71 * ln 2)
# => D/K must approximate log2(3) with error < 1/(3*2^71*ln2) ~ 4e-22 / 1  per unit K
lg = math.log2(3)
convs = []
frac = lg
a0 = int(frac); h0, k0 = a0, 1; h1, k1 = 1, 0
x = frac
cf = []
for _ in range(30):
    ai = int(x)
    cf.append(ai)
    x = 1 / (x - ai) if x != ai else 0
    h0, h1 = ai * h0 + h1, h0
    k0, k1 = ai * k0 + k1, k0
    convs.append((h0, k0, abs(lg - h0 / k0)))
    if x == 0: break
# minimum K such that |D/K - lg| can be < eps = 1/(3*2^71*ln2*K):
eps_per = 1 / (3 * (2 ** 71) * math.log(2))
kmin = None
for h, kq, err in convs:
    # best possible error at denominator kq is ~1/(kq*k_{q+1}); need err*kq^... simpler: err < eps_per*...
    if err < eps_per:
        kmin = kq
        break
out["log2_3_convergents"] = [(h, kq, f"{err:.3e}") for h, kq, err in convs[:16]]
out["needed_precision_per_K"] = f"{eps_per:.3e}"
out["min_odd_steps_K_from_convergents"] = kmin

# ============ E7: needle DP to t=120 + survivor structure ============
alpha = 1 / math.log2(3)
surv_frac = {}
dp = {0: 1}
for j in range(1, 121):
    nd = {}
    for K, c in dp.items():
        for b in (0, 1):
            K2 = K + b
            if K2 >= alpha * j - 1e-12:
                nd[K2] = nd.get(K2, 0) + c
    dp = nd
    if j in (20, 40, 60, 80, 100, 120):
        surv_frac[j] = sum(dp.values()) / 2 ** j
out["needle_survivor_fraction"] = {j: f"{v:.3e}" for j, v in surv_frac.items()}
kl = alpha * math.log2(2 * alpha) + (1 - alpha) * math.log2(2 * (1 - alpha))
out["thinning_rate_bits_per_step"] = round(kl, 4)

# survivor structure: actual integers < 2^24 that don't drop below self in 60 T-steps
def survives(n, t):
    m = n
    for _ in range(t):
        m = (3 * m + 1) // 2 if m % 2 else m // 2
        if m < n: return False
    return True
survivors = [n for n in range(3, 1 << 22, 2) if survives(n, 40)]
out["n_survivors_lt_2^22_40steps"] = len(survivors)
mean_t1_surv = sum(trailing_ones(n) for n in survivors) / len(survivors)
mean_t1_all = sum(trailing_ones(n) for n in range(3, 1 << 22, 2)) / len(range(3, 1 << 22, 2))
out["mean_trailing1s"] = dict(survivors=round(mean_t1_surv, 3), all_odd=round(mean_t1_all, 3))

# figure: needle
xs = sorted(surv_frac)
plt.figure(figsize=(7, 4.2))
plt.semilogy(xs, [surv_frac[j] for j in xs], "o-")
plt.xlabel("T-steps t"); plt.ylabel("fraction of residues mod 2^t still alive")
plt.title(f"E7: never-drop needle thins at 2^(-{kl:.3f} t) — vanishing, never empty")
plt.tight_layout()
plt.savefig(os.path.join(FIG, "e7_needle.png"), dpi=110)

print(json.dumps(out, indent=1))
