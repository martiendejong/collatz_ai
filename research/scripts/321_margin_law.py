"""
Pattern: the certified margins R - c2/c0 across (lam, k).
New exact identity (via mu2 = R*mu0): margin = (g2 - R*g0)/c0.
Prediction: margin decays in k at rate sqrt(c(lam)) (the amplitude rate),
since g-scale ~ relgap ~ sqrt(c)^k and (G-1) stabilizes.
Data: 72 rigorous margins from cert_lemmaA_k5-10.txt.
"""
import re
import numpy as np

txt = open('E:/projects/collatz/research/certificates/cert_lemmaA_k5-10.txt').read()
pat = re.compile(r"lam=(\S+) k=(\d+): rho=\[([0-9.]+),.*?\n.*?R-c2/c0=\[\+?([0-9.e-]+)")
data = {}
for m in pat.finditer(txt):
    lam_s, k, rho, marg = m.group(1), int(m.group(2)), float(m.group(3)), float(m.group(4))
    data.setdefault(lam_s, {})[k] = (marg, rho)

print("lam | margins k=5..10 | ratio(k)/(k-1) | mean rate | sqrt(c) pred")
# c(lam) estimates from today's instruments: c ~ measured deep prefix-plateau
c_est = {"21/20": 0.41, "13/10": 0.55, "17/10": 0.70, "2": 0.824}  # 13/10 rough
for lam_s in ["21/20", "11/10", "6/5", "13/10", "7/5", "3/2", "8/5", "17/10", "9/5", "19/10", "39/20", "2"]:
    if lam_s not in data: continue
    d = data[lam_s]
    ks = sorted(d)
    margs = [d[k][0] for k in ks]
    rats = [margs[i+1]/margs[i] for i in range(len(margs)-1)]
    pred = f"{np.sqrt(c_est[lam_s]):.3f}" if lam_s in c_est else "  -  "
    print(f"{lam_s:>6} | " + " ".join(f"{m:.2e}" for m in margs) + " | " +
          " ".join(f"{r:.3f}" for r in rats) + f" | {np.mean(rats[-3:]):.3f} | {pred}")

# verify the margin identity on cached vectors
print("\nExact identity check: margin = (g2 - R*g0)/c0")
CACHE = "E:/projects/collatz/research/cache"
from math import log2
for lam, k in [(1.05, 13), (1.70, 13)]:
    v = np.load(f"{CACHE}/v_lam{lam:.2f}_k{k}.npy")
    rho = float(open(f"{CACHE}/rho_lam{lam:.2f}_k{k}.txt").read())
    t = lam**-2/rho
    R = (t*t+lam)/(1+t*lam)
    N = v.size; Nl = N//3
    cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
    mu0 = v[0::3].mean(); mu2 = v[2::3].mean()
    c0 = cb[0::3].mean(); c2 = cb[2::3].mean()
    g0 = mu0-c0; g2 = mu2-c2
    lhs = R - c2/c0
    rhs = (g2 - R*g0)/c0
    print(f"lam={lam} k={k}: R-c2/c0={lhs:.8e}  (g2-R*g0)/c0={rhs:.8e}  ratio={lhs/rhs:.10f}")
