import numpy as np
from math import log2

L = log2(3.0)  # 1.58496...

# measured sequences at lambda=1.05 (Obs 484/485/488/491, Scripts 281/282/285/289)
ks = np.arange(3, 21)
ratio = np.array([3.62408, 1.66387, 1.23605, 1.20471, 1.11916, 1.15504, 1.10348,
                  1.09917, 1.08878, 1.08945, 1.08554, 1.08380, 1.08300, 1.08217,
                  1.07867, 1.07856, 1.07794, 1.07769])
# s2/s0 at lambda=1.05 (Scripts 291/292), k=4..17
ks_s = np.arange(4, 18)
s2s0 = np.array([0.28711, 0.45195, 0.58403, 0.74271, 0.55419, 0.66033, 0.67704,
                 0.70178, 0.70344, 0.71629, 0.72299, 0.72815, 0.72887, 0.74094])

def resonance(k):
    x = k * L
    e = x - round(x)          # signed: >0 means 3^k just above 2^m
    return e

print("k | k*log2(3) | signed e(k) | |e| | nearest 2^m | ratio-anomaly (2nd diff) | s2/s0 anomaly")
for k in range(3, 22):
    e = resonance(k)
    m = round(k*L)
    a_r = ""
    if 4 <= k <= 19:
        i = k - 3
        a = ratio[i+1] - 2*ratio[i] + ratio[i-1]
        a_r = f"{a:+.5f}"
    a_s = ""
    if 5 <= k <= 16:
        j = k - 4
        a2 = s2s0[j+1] - 2*s2s0[j] + s2s0[j-1]
        a_s = f"{a2:+.5f}"
    flag = " <== RESONANT" if abs(e) < 0.10 else ""
    print(f"{k:2d} | {k*L:8.4f} | {e:+.3f} | {abs(e):.3f} | 2^{m} | {a_r:>9s} | {a_s:>9s}{flag}")

# correlation in the smooth regime k=9..19 for ratio second-diffs
kk = np.arange(9, 20)
second = np.array([ratio[i+1] - 2*ratio[i] + ratio[i-1] for i in kk-3])
absres = np.array([abs(resonance(int(k))) for k in kk])
sgnres = np.array([resonance(int(k)) for k in kk])
from scipy import stats as st
try:
    r1 = st.spearmanr(absres, np.abs(second))
    r2 = st.spearmanr(sgnres, second)
    print(f"\nSpearman |e| vs |2nd-diff ratio| (k=9..19): rho={r1.statistic:+.3f} p={r1.pvalue:.3f}")
    print(f"Spearman signed e vs signed 2nd-diff:      rho={r2.statistic:+.3f} p={r2.pvalue:.3f}")
except Exception as ex:
    # fallback without scipy
    def spearman(a, b):
        ra = np.argsort(np.argsort(a)); rb = np.argsort(np.argsort(b))
        return np.corrcoef(ra, rb)[0,1]
    print(f"\nSpearman |e| vs |2nd diff|: {spearman(absres, np.abs(second)):+.3f}")
    print(f"Spearman signed:            {spearman(sgnres, second):+.3f}")

# s2/s0 correlations k=6..16
kk2 = np.arange(6, 17)
second2 = np.array([s2s0[j+1] - 2*s2s0[j] + s2s0[j-1] for j in kk2-4])
absres2 = np.array([abs(resonance(int(k))) for k in kk2])
try:
    r3 = st.spearmanr(absres2, np.abs(second2))
    print(f"Spearman |e| vs |2nd-diff s2/s0| (k=6..16): rho={r3.statistic:+.3f} p={r3.pvalue:.3f}")
except Exception:
    print(f"Spearman s2/s0: {spearman(absres2, np.abs(second2)):+.3f}")
