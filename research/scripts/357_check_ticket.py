# 357: checking the 186-billion ticket (K=72,057,431,991, N=114,208,327,604)
# as far as is feasible. Full check = deciding (2^N - 3^K) | c over C(N-1,K-1)
# compositions: world-record territory, not computable. Feasible checks:
#  (1) exact window: Delta = N ln2 - K ln3, required sum 1/(3 n_i) = Delta,
#      kill threshold V* = K/(3 Delta) (harmonic mean of elements).
#      Current verification V = 2^71 (Barina, Jan 2025): alive iff V < V*.
#  (2) phase forcing (new for this ticket): per climb from phase-minimum n0 the
#      contributions sum to <= C/n0 (geometric); with all minima > V this forces
#      m >= Delta * V / C phases. Compare with Hercher's m >= 92.
#  (3) the death schedule: how many extra verification bits kill the ticket,
#      and when at the project's historical pace.
from decimal import Decimal, getcontext
from math import log2
getcontext().prec = 80

K = 72057431991
N = 114208327604
ln2 = Decimal(2).ln()
ln3 = Decimal(3).ln()
Delta = N * ln2 - K * ln3          # = sum ln(1+1/(3 n_i)) ~ sum 1/(3 n_i)
Vstar = K / (3 * Delta)            # required harmonic mean of elements
V = Decimal(2) ** 71               # current verification bound

print(f"ticket: K = {K:,} odd steps, N = {N:,} halvings, length {K+N:,}")
print(f"exact window Delta = N ln2 - K ln3 = {Delta:.6E}")
print(f"required sum of 1/(3 n_i)          = {Delta:.6E}  (all pushes positive)")
print(f"kill threshold V* (harmonic mean)  = {Vstar:.6E} = 2^{float(Vstar.ln()/ln2):.4f}")
print(f"current verification V = 2^71      = {V:.6E}")
alive = V < Vstar
print(f"STATUS: {'ALIVE (V < V*)' if alive else 'DEAD (V >= V*)'}, slack factor {float(Vstar/V):.3f}")
print()
# (2) phase forcing
for Cc, label in [(Decimal(1), "climb-only assignment (C=1)"),
                  (Decimal(7) / 3, "climb+descent assignment (C=7/3, conservative)")]:
    m_min = Delta * V / Cc
    print(f"phase lower bound, {label}: m >= {float(m_min):.3E}")
    print(f"   -> average phase length <= {float(K / m_min):.1f} odd steps "
          f"(Hercher 2022 requires m >= 92: satisfied but vastly superseded)")
print()
# (3) death schedule
import math
bits_needed = float((Vstar / V).ln() / ln2)
# project pace: 2^68 (2021) -> 2^71 (Jan 2025): ~0.79 bits/yr
pace = (71 - 68) / 3.8
print(f"extra verification bits needed to kill the ticket: {bits_needed:.3f}")
print(f"historical pace ~{pace:.2f} bits/yr -> expected death in ~{bits_needed/pace:.1f} yr from Jan 2025")
print()
# sanity: the next ticket after death
K2, N2 = 137528045312, 217976794617
D2 = N2 * ln2 - K2 * ln3
V2 = K2 / (3 * D2)
print(f"next ticket after death: K = {K2:,}, kill threshold 2^{float(V2.ln()/ln2):.2f}")
