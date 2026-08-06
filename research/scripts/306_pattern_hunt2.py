import numpy as np
from math import log2

# ============ A: gamma-ladder decay structure ============
ks = np.array([9, 11, 12, 13, 15, 17, 19, 20, 21])
gam = np.array([0.8168, 0.8418, 0.852, 0.8624, 0.8805, 0.8953, 0.9069, 0.9146, 0.9184])
one_m = 1 - gam
print("=== A: 1-gamma geometric decay ===")
for i in range(1, len(ks)):
    dk = ks[i]-ks[i-1]
    r = (one_m[i]/one_m[i-1])**(1/dk)
    print(f"k={ks[i-1]}->{ks[i]}: rate/step = {r:.4f}")
lr = np.polyfit(ks, np.log(one_m), 1)
rate = np.exp(lr[0])
resid = np.log(one_m) - np.polyval(lr, ks)
print(f"global fit: (1-gamma) ~ C * {rate:.4f}^k, max |log-resid| = {np.abs(resid).max():.4f}")
k95 = 21 + np.log(0.05/one_m[-1])/np.log(rate)
print(f"extrapolated gamma=0.95 crossing: k ~ {k95:.1f} (fork: DENSITY predicts ~27+-3, CEILING never)")

# ============ B: rho convergence rate vs cascade ratio ============
print("\n=== B: rho(k) convergence at lam=1.05 vs sqrt(cascade plateau) ===")
rhos = {8: 1.573725431, 10: 1.575493204, 12: 1.576193953, 13: 1.576366173, 14: 1.576475227}
r1213 = rhos[13]-rhos[12]; r1314 = rhos[14]-rhos[13]
r1012 = (rhos[12]-rhos[10])/2; r0810 = (rhos[10]-rhos[8])/2
print(f"increments: 8-10:{r0810:.3e}/step 10-12:{r1012:.3e}/step 12-13:{r1213:.3e} 13-14:{r1314:.3e}")
print(f"ratios: {r1012/r0810:.3f} (2-step avg) {r1213/r1012:.3f} {r1314/r1213:.3f}")
print(f"sqrt(prefix plateau 0.41) = {np.sqrt(0.41):.3f}  <- match?")

# ============ C: digit-1 profile predicted from exact identities ============
print("\n=== C: digit-1 main effect vs analytical prediction ===")
ALPHA = log2(3.0)
for lam in [1.05, 1.70]:
    k = 13
    v = np.load(f"E:/projects/collatz/research/cache/v_lam{lam:.2f}_k{k}.npy")
    rho = float(open(f"E:/projects/collatz/research/cache/rho_lam{lam:.2f}_k{k}.txt").read())
    A = lam**-2; B1 = lam**(ALPHA-2); B3 = lam**(ALPHA-1)
    N = v.size; Nl = N//3
    cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
    F = np.log2(v); F -= F.mean()
    i = np.arange(N)
    d1 = (i//3) % 3
    # measured digit-1 profile (overall and per class r = i mod 3)
    meas = np.array([F[d1==d].mean() for d in range(3)])
    meas -= meas.mean()
    # prediction: class-r rows receive cb sub-class q(r, m mod 3):
    # r=0: q = (4m) mod 3 = m mod 3 = d1 ; r=2: q = (2m+1) mod 3 = transposition(0<->1) of d1
    # log-shift ~ (B_r * (c_q - cbar)) / (rho * v-scale) -> profile proportional to c-spread
    cq = np.array([cb[q::3].mean() for q in range(3)])
    cbar = cq.mean()
    mu = v.mean()
    # per-class contributions to the OVERALL digit-1 effect (log-linearized):
    pred = np.zeros(3)
    for d in range(3):
        # class 0 rows with m mod 3 = d get c_d ; class 2 rows get c_{tau(d)}, tau = (0<->1)
        tau = [1,0,2][d]
        contrib0 = B1*(cq[d]-cbar)/(rho*mu*np.log(2)*3)   # 1/3 of rows are class 0
        contrib2 = B3*(cq[tau]-cbar)/(rho*mu*np.log(2)*3)
        pred[d] = contrib0 + contrib2
    pred -= pred.mean()
    cos = pred@meas/np.linalg.norm(pred)/np.linalg.norm(meas)
    scale = np.linalg.norm(meas)/np.linalg.norm(pred)
    print(f"lam={lam}: measured d1-profile {meas.round(5)} predicted-shape {pred.round(6)}")
    print(f"   cosine(pred, meas) = {cos:+.4f}  amplitude ratio meas/pred = {scale:.3f}")
