"""Pattern: digit-profiles of the cb-field vs the v-field (on [0,Nl)).
If parallel (cos ~ 1) with amplitude ratio a_q, the c-hierarchy is v-slaved
and a_q is a new fixed profile; a_q relates to the min-transfer c_q of Obs 504."""
import numpy as np
from math import log2
CACHE = "E:/projects/collatz/research/cache"

for lam in [1.05, 1.70]:
    k = 13
    v = np.load(f"{CACHE}/v_lam{lam:.2f}_k{k}.npy")
    N = v.size; Nl = N//3
    cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
    vbar = (v[:Nl]+v[Nl:2*Nl]+v[2*Nl:])/3.0
    j = np.arange(Nl, dtype=np.int64)
    P = k-2
    digs = np.empty((P, Nl), dtype=np.int8)
    x = j.copy()
    for p in range(P):
        digs[p] = x % 3; x //= 3
    print(f"\n=== lam={lam} k={k} ===")
    print("p | cos(cb_p, vbar_p) | amp ratio cb/vbar | amp ratio (1 - ...) vs c_q")
    for p in range(min(P, 8)):
        prof_cb = np.array([cb[digs[p]==d].mean() for d in range(3)]); prof_cb -= prof_cb.mean()
        prof_vb = np.array([vbar[digs[p]==d].mean() for d in range(3)]); prof_vb -= prof_vb.mean()
        ncb, nvb = np.linalg.norm(prof_cb), np.linalg.norm(prof_vb)
        cos = prof_cb@prof_vb/(ncb*nvb) if ncb>0 and nvb>0 else 0
        print(f"{p} | {cos:+.4f} | {ncb/nvb:.4f}")
