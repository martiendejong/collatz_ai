"""R256: SPATIAL CORRELATION STRUCTURE of the K-L eigenprofile at k=13/15.
The lattice coefficients (a,c) emerge from spatial correlations (R206-210);
this measures that structure directly: autocorrelation of the triple-CV field
over 3-adic class index, its correlation length, and CV-vs-mean coupling.
Triple(i) = C[i], C[i+3^{k-2}], C[i+2*3^{k-2}] (children of parent class 3i+2)."""
import sys
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")

for k, path in ((13, "certificates/cert_k13.npy"), (15, "certificates/cert_k15.npy")):
    C = np.load(path, mmap_mode="r").astype(np.float64)
    P = 3 ** (k - 2)
    T = np.stack([C[:P], C[P:2*P], C[2*P:3*P]])
    mean = T.mean(0)
    cv = T.std(0) / mean
    print(f"k={k}: {P:,} parents  CV mean {cv.mean():.4f}  std {cv.std():.4f}")
    x = cv - cv.mean()
    v = (x*x).mean()
    print("  lag autocorrelation (3-adic index):")
    row = []
    for lag in (1, 2, 3, 9, 27, 81, 243, 729, 2187, 6561, 19683, 59049):
        r = (x[:-lag]*x[lag:]).mean()/v
        row.append(f"lag {lag}: {r:+.4f}")
    print("   " + " | ".join(row[:6]))
    print("   " + " | ".join(row[6:]))
    # correlation length: first lag where |r| < 1/e
    for lag in range(1, 5000):
        if abs((x[:-lag]*x[lag:]).mean()/v) < 1/np.e:
            print(f"  correlation length xi = {lag} (first |r| < 1/e)")
            break
    mm = mean - mean.mean()
    r_cm = (x*mm).mean()/np.sqrt(v*(mm*mm).mean())
    print(f"  corr(CV, triple-mean) = {r_cm:+.4f}")
    # is the CV field self-similar? block-average at scale 3 and re-measure lag-1
    cb = cv[:3*(P//3)].reshape(-1, 3).mean(1)
    xb = cb - cb.mean()
    r1b = (xb[:-1]*xb[1:]).mean()/(xb*xb).mean()
    r1 = (x[:-1]*x[1:]).mean()/v
    print(f"  lag-1 raw {r1:+.4f} vs after 3-blocking {r1b:+.4f} (self-similarity test)")
    print()
