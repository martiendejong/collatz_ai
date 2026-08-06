"""
Induction structure for Lemma alpha at deeper q.
1. EXACT fibered class-1 identity: E[v1|cell] = t * E[v0|mapped cell] per cell, every level.
2. Per-cell local structure: within each level-q cell, the 3 sub-cells give a local
   mean-triple u_C and gap-triple w_C. Chebyshev per cell needs similar ordering.
   Measure: fraction of cells with positive local cov, aggregate decomposition,
   and self-similarity of the average local shape.
"""
import numpy as np
from math import log2

ALPHA = log2(3.0)
CACHE = "E:/projects/collatz/research/cache"

for lam in [1.05, 1.70]:
    k = 13
    v = np.load(f"{CACHE}/v_lam{lam:.2f}_k{k}.npy")
    rho = float(open(f"{CACHE}/rho_lam{lam:.2f}_k{k}.txt").read())
    A = lam**-2
    t = A/rho
    N = v.size; Nl = N//3
    v0 = v[0::3]; v1 = v[1::3]
    # 1. fibered class-1 identity per cell at level q: E[v1 | s in cell] vs t*E[v0 | mapped]
    # v1[s] = t*v0[(4s+2) mod Nl] exactly pointwise (already known) -> cell version trivial.
    # The nontrivial fibered check: does the local sub-triple of v0-cells map to v1 sub-triples
    # with the SAME digit alignment? (4s+2 mod 3^q is a bijection on cells: yes, affine.)
    s = np.arange(Nl)
    err = np.abs(v1 - t*v0[(4*s+2) % Nl]).max()
    print(f"\n=== lam={lam} ===  pointwise class-1 identity err {err:.1e} (fibered version follows)")

    # 2. per-cell local (u, w) structure on the cb index space
    cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
    vbar3 = (v[:Nl] + v[Nl:2*Nl] + v[2*Nl:])/3.0
    G = vbar3 - cb
    P = k-2
    for q in [1, 2, 4, 6]:
        M = 3**(q+1)          # sub-cell count at split level
        Mp = 3**q             # parent cells
        # conditional means on sub-cells (level q) as (parent, subdigit) arrays
        cmV = vbar3.reshape(Nl//M, M).mean(axis=0)   # index = low q+1 digits
        cmG = G.reshape(Nl//M, M).mean(axis=0)
        # reshape: sub-cell index = parent + d*3^q  (digit q is the HIGH digit of the low-block)
        U = cmV.reshape(3, Mp).T   # [parent, d]
        W = cmG.reshape(3, Mp).T
        Uc = U - U.mean(axis=1, keepdims=True)
        Wc = W - W.mean(axis=1, keepdims=True)
        loccov = (Uc*Wc).mean(axis=1)
        frac_pos = (loccov > 0).mean()
        agg = loccov.mean()
        # weight: parent-level b comes from aggregate of local cov + between-parent part
        # self-similar shape: average normalized local u-triple
        nrm = np.linalg.norm(Uc, axis=1, keepdims=True) + 1e-300
        avg_shape_u = (Uc/nrm).mean(axis=0)
        nrmw = np.linalg.norm(Wc, axis=1, keepdims=True) + 1e-300
        avg_shape_w = (Wc/nrmw).mean(axis=0)
        # coherence: how aligned are local shapes across cells (1 = identical shape)
        coh_u = np.linalg.norm(avg_shape_u)
        coh_w = np.linalg.norm(avg_shape_w)
        print(f"q={q}: frac cells loccov>0: {frac_pos:.3f} | agg loccov {agg:+.2e} | "
              f"shape coherence u={coh_u:.3f} w={coh_w:.3f} | avg u-shape {avg_shape_u.round(3)} w-shape {avg_shape_w.round(3)}")
