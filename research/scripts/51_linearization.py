"""R409-420: THE LINEARIZATION LOOP, step 1: the min-mean gap law.
For triples (x1,x2,x3): min = mean*(1 - c1*CV + c2*CV^2 + ...)?
For Gaussian triples c1 = 3/(2*sqrt(pi))? -- measure empirically on the actual
certificate triples at every level p and every k. If c1 is stable and the
quadratic term small, the K-L min-operator linearizes as
   (K-L update) = (roulette-weighted mean) * (1 - c1*CV_local)
and the CV-cascade contraction (theta<1) closes the alpha->1 loop."""
import sys, math
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")

for k, path in ((13, "certificates/cert_k13.npy"), (17, "certificates/cert_k17.npy")):
    C = np.load(path, mmap_mode="r")
    N = 3 ** (k - 1)
    print(f"\n=== k={k} ===")
    print(f"{'level p':>7} {'mean CV':>8} {'c1 (fit)':>9} {'c2 (fit)':>9} {'R2':>7}")
    for p in (0, 2, 4, 6, 8):
        B = 3 ** p
        CH = 3 ** 13
        cvs, gaps = [], []
        for lo in range(0, min(N, 3*CH*(3**p)), 3 * B * (CH // max(1, B // 3**6))):
            idx = np.arange(lo, min(lo + 3*B*200, N), dtype=np.int64)
            sel = idx[(idx // B) % 3 == 0]
            sel = sel[sel + 2 * B < N][:200000]
            if sel.size == 0: continue
            t = np.stack([C[sel], C[sel + B], C[sel + 2 * B]]).astype(np.float64)
            mn = t.mean(0)
            cv = t.std(0) / mn
            gap = 1 - t.min(0) / mn
            cvs.append(cv); gaps.append(gap)
            if sum(x.size for x in cvs) > 300000: break
        cv = np.concatenate(cvs); gap = np.concatenate(gaps)
        ok = cv > 1e-9
        cv, gap = cv[ok], gap[ok]
        A = np.column_stack([cv, cv**2])
        coef, *_ = np.linalg.lstsq(A, gap, rcond=None)
        pred = A @ coef
        r2 = 1 - ((gap - pred)**2).sum()/((gap - gap.mean())**2).sum()
        print(f"{p:>7} {cv.mean():>8.4f} {coef[0]:>9.4f} {coef[1]:>9.4f} {r2:>7.4f}")
print("\n[Gaussian iid triple reference: E[1-min/mean] = c1*CV with c1 = 1.5*E|z|-ish ~ 1.1958/",
      f"{3/(2*math.sqrt(math.pi)):.4f} = 3/(2 sqrt pi)]")
