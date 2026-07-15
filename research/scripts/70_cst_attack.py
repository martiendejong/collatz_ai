"""R711-725: TERRAS CST VERIFICATION + MARGINS. CST: tau(n) (first t with
coefficient 3^u < 2^t) equals sigma(n) (first actual drop) for all n.
A violation needs value a*n+b >= n at t=tau, i.e. n < b/(1-a) with n in the
class r mod 2^tau. Enumerate all tau=t classes (t <= 24), compute threshold
b/(1-a), compare with minimal class member r: violation possible only if
r < threshold -> check those n directly. Also: margin statistics."""
import sys, math
from fractions import Fraction
sys.stdout.reconfigure(encoding="utf-8")

# extend survivor-style: states (r, u, j, b) with 3^u >= 2^j so far (not yet stopped)
cur = [(1, 0, 0, 0)]
worst_margin = []
viol = 0
checked = 0
for t in range(1, 25):
    nxt = []
    M = 1 << t
    stopped_here = []
    for r, u, j, b in cur:
        for r2 in (r, r + (1 << (t - 1))):
            val = (pow(3, u, M << 2) * r2 + b) >> j
            if val & 1:
                u2, b2 = u + 1, 3 * b + (1 << j)
            else:
                u2, b2 = u, b
            j2 = j + 1
            if 3 ** u2 < (1 << j2):
                stopped_here.append((r2, u2, j2, b2))
            else:
                nxt.append((r2, u2, j2, b2))
    cur = nxt
    for r2, u2, j2, b2 in stopped_here:
        # tau = j2 for all n == r2 mod 2^j2; threshold n < b2/(2^j2 - 3^u2) * ... :
        # value after j2 steps = (3^u2 n + b2)/2^j2 >= n  <=>  n <= b2/(2^j2-3^u2)
        thr = Fraction(b2, (1 << j2) - 3 ** u2)
        checked += 1
        if r2 <= thr:
            # potential violation: check n = r2 directly
            n = r2; m = n; drop = False
            for _ in range(j2):
                m = 3*m+1 if m % 2 else m // 2
                if m < n: drop = True; break
            if not drop: viol += 1; print(f"  VIOLATION: n={r2}, tau={j2}")
        worst_margin.append(float(Fraction(r2) / thr) if thr > 0 else 1e9)
print(f"t <= 24: {checked:,} stopping classes checked, violations = {viol}")
wm = sorted(worst_margin)[:8]
print("smallest margins r/threshold (must stay > 1):", [f"{x:.2f}" for x in wm])
