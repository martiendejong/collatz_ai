# 389: THE COUPLING OF THE TWO CLOCKWORKS — the first direct measurement of
# the self-referential sampling. Head = archimedean angle phi = frac(log2 n)
# (the Weyl rotation coordinate); tail = 2-adic clock theta (discrete 3-log).
# Under no-conspiracy the joint distribution should be uniform-independent.
# Any head-tail correlation would be the first measurable trace of the
# coupling on which the wall lives.
import numpy as np
import random
from math import log2

M = 16
MOD = 1 << M
ORD = 1 << (M-2)
tab = {}
x = 1
for j in range(ORD):
    tab[x] = j; x = (x*3) % MOD
def theta(n):
    r = n % MOD
    return tab[r] if r in tab else tab[(-r) % MOD]

def head_angle(n):
    L = n.bit_length()
    top = n >> (L - 53) if L > 53 else n
    return log2(top) - (top.bit_length() - 1) if False else (log2(top) + (L - top.bit_length())) % 1.0

random.seed(389)
NB = 16
H = np.zeros((NB, NB))
tot = 0
mi_pairs = []
for _ in range(400):
    n = random.getrandbits(700) | 1
    for _ in range(1200):
        m_ = 3*n + 1
        v = (m_ & -m_).bit_length() - 1
        n = m_ >> v
        if n.bit_length() < 80:
            break
        ph = head_angle(n)
        th = theta(n)
        H[int(ph*NB) % NB, (th >> (M-2-4)) % NB] += 1
        tot += 1
P = H/tot
Pr = P.sum(axis=1, keepdims=True); Pc = P.sum(axis=0, keepdims=True)
chi2 = float(((H - tot*Pr*Pc)**2/(tot*Pr*Pc + 1e-12)).sum())
df = (NB-1)**2
mi = float((P*np.log2(P/(Pr*Pc) + 1e-300)).sum())
# uniformity of each marginal
zr = float(np.abs(H.sum(axis=1) - tot/NB).max()/np.sqrt(tot/NB))
zc = float(np.abs(H.sum(axis=0) - tot/NB).max()/np.sqrt(tot/NB))
print(f"{tot} stappen; 16x16 kop-hoek x staart-klok")
print(f"marginalen uniform: kop max|dev| {zr:.2f} sd, staart {zc:.2f} sd")
print(f"onafhankelijkheid: chi2 = {chi2:.0f} bij {df} df (verwacht {df} +- {int((2*df)**0.5)})")
print(f"wederzijdse informatie: {mi:.2e} bits (0 = onafhankelijk)")
