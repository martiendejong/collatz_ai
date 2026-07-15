"""R577-585: PROFILE COLLAPSE TEST. CV_p(k) for all levels p, k=13/15/17
(k=20 already computed). Two alignments: from the fine end (p) and from the
top (d = k-2-p). If profiles collapse in one alignment -> traveling wave;
the alpha->1 question reduces to the wave's amplitude law."""
import sys
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")

prof20 = [0.4629,0.3691,0.3145,0.2754,0.2402,0.2092,0.1806,0.1553,0.1333,
          0.1145,0.0978,0.0828,0.0695,0.0578,0.0473,0.0371,0.0277,0.0195]
profs = {20: prof20}
for k, path in ((13,"certificates/cert_k13.npy"),(15,"certificates/cert_k15.npy"),
                (17,"certificates/cert_k17.npy")):
    C = np.load(path, mmap_mode="r")
    N = 3 ** (k - 1)
    CH = 3 ** 13
    out = []
    for p in range(1, k - 1):
        B = 3 ** p
        acc = 0.0; cnt = 0
        for lo in range(0, N, 3 * CH):
            hi = min(lo + 3 * CH, N)
            idx = np.arange(lo, hi, dtype=np.int64)
            sel = idx[(idx // B) % 3 == 0]
            sel = sel[sel + 2 * B < N]
            if sel.size == 0: continue
            t = np.stack([C[sel], C[sel + B], C[sel + 2 * B]]).astype(np.float64)
            cv = t.std(0) / t.mean(0)
            acc += cv.sum(); cnt += cv.size
        out.append(acc / cnt)
    profs[k] = out
    print(f"k={k}: " + " ".join(f"{x:.4f}" for x in out))

print("\nfine-end alignment CV_p(k) (rows p=1..6):")
for p in range(6):
    print(f"  p={p+1}: " + " ".join(f"k{k}:{profs[k][p]:.4f}" for k in (13,15,17,20)))
print("\ntop alignment CV at d=k-2-p (rows d=0..5, d=0 = top):")
for d in range(6):
    row = []
    for k in (13,15,17,20):
        arr = profs[k]
        row.append(f"k{k}:{arr[len(arr)-1-d]:.4f}")
    print(f"  d={d}: " + " ".join(row))
print("\nper-level decay ratio (k=20, fine->coarse):",
      " ".join(f"{prof20[i+1]/prof20[i]:.3f}" for i in range(0, 17, 2)))
