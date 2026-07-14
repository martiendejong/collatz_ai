"""E26: period-5 realization census, k,w <= 5 (9.77M shapes, sign-free)."""
import sys, json
sys.stdout.reconfigure(encoding="utf-8")

def compose5(shape):
    A, B, S = 1, 0, 0
    for i in range(5):
        k, w = shape[i]
        kn = shape[(i + 1) % 5][0]
        A *= 3 ** k
        B = B * 3 ** k + ((1 << w) - 1) * (1 << S)
        S += w + kn
    return A, B, S

found = []
R = range(1, 6)
p3 = [3 ** i for i in range(6)]
for k1 in R:
 for w1 in R:
  for k2 in R:
   for w2 in R:
    for k3 in R:
     for w3 in R:
      for k4 in R:
       for w4 in R:
        for k5 in R:
         for w5 in R:
            shape = ((k1, w1), (k2, w2), (k3, w3), (k4, w4), (k5, w5))
            A, B, S = compose5(shape)
            D = (1 << S) - A
            if D == 0: continue
            if B % abs(D) == 0:
                a = B // D
                if a % 2 != 0 and a != 0:
                    if len(set(shape)) > 1 or a not in (1, -1):
                        found.append((shape, a))
print(json.dumps(dict(checked=5 ** 10,
    nontrivial=[f for f in found if not all(x == (1, 1) for x in f[0])][:10],
    total_hits=len(found)), indent=1))
