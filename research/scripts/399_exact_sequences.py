# 399: pattern recognition on EXACT integer sequences (cheap, memory-light).
# (A) Mirror sequence: A_j = N_j / 2^j with N_j integer (exact enumeration,
#     j <= 20 to stay memory-safe next to the background job). Search a
#     linear recurrence with rational coefficients on A_j, and separately on
#     the even/odd subsequences (period-2 structure).
# (B) Bias channel: exact conditional biases of f_t on the best motif slice,
#     full enumeration mod 2^(t+1) (tiny) — exact dyadics, sequence identify.
import numpy as np
from fractions import Fraction

# ---- (A) exact mirror numerators ----
As = {}
for j in range(1, 21):
    Nj = 1 << (j+1)
    n = np.arange(Nj, dtype=np.int64)
    for _ in range(j):
        odd = (n & 1).astype(bool)
        n = np.where(odd, (3*n + 1) >> 1, n >> 1)
    s = 1 - 2*(n & 1).astype(np.int64)
    idx = np.arange(Nj)
    neg = (Nj - idx) % Nj
    m = (idx & 1).astype(bool)
    tot = int((s[m]*s[neg[m]]).sum())
    As[j] = Fraction(tot, 1 << j)
    del n, s, idx, neg, m
print("exacte A_j (teller/2^j):")
for j in sorted(As):
    f = As[j]
    print(f"  j={j:2d}: {f}  = {float(f):+.6f}")

def find_recurrence(seq, max_order=6):
    # seq: list of Fractions; find minimal r with a_n = sum c_i a_{n-i}
    for r in range(1, max_order+1):
        if len(seq) < 2*r + 2:
            break
        # solve from r equations
        import itertools
        rows = [[seq[k+i] for i in range(r)] for k in range(r)]
        rhs = [seq[k+r] for k in range(r)]
        # Gaussian elimination over fractions
        A = [row[:] + [b] for row, b in zip(rows, rhs)]
        nr = len(A)
        ok = True
        for col in range(r):
            piv = next((rw for rw in range(col, nr) if A[rw][col] != 0), None)
            if piv is None:
                ok = False; break
            A[col], A[piv] = A[piv], A[col]
            inv = 1/A[col][col]
            A[col] = [x*inv for x in A[col]]
            for rw in range(nr):
                if rw != col and A[rw][col] != 0:
                    f2 = A[rw][col]
                    A[rw] = [x - f2*y for x, y in zip(A[rw], A[col])]
        if not ok:
            continue
        cs = [A[i][r] for i in range(r)]
        # verify on the rest
        good = all(sum(cs[i]*seq[k+i] for i in range(r)) == seq[k+r]
                   for k in range(len(seq)-r))
        if good:
            return r, cs
    return None, None

seq_all = [As[j] for j in sorted(As)]
seq_even = [As[j] for j in sorted(As) if j % 2 == 0]
seq_odd = [As[j] for j in sorted(As) if j % 2 == 1]
for name, seq in [("volledig", seq_all), ("even j", seq_even), ("oneven j", seq_odd)]:
    r, cs = find_recurrence(seq)
    if r:
        print(f"recurrentie ({name}): orde {r}, coeff {[str(c) for c in cs]}")
    else:
        print(f"recurrentie ({name}): geen t/m orde 6")

# ---- (B) exact bias channel dyadics ----
print("\nexacte bias-kanaal-tabel (volledige enumeratie):")
print(f"{'t':>3} {'bias(f_t)':>12} {'beste (a,b)':>11} {'bias|ab=1 exact':>16} {'bias|ab=0 exact':>16}")
for t in range(4, 15):
    Nt = 1 << (t+1)
    n = np.arange(Nt, dtype=np.int64)
    for _ in range(t):
        odd = (n & 1).astype(bool)
        n = np.where(odd, (3*n + 1) >> 1, n >> 1)
    p = (n & 1).astype(np.int64)
    x = np.arange(Nt)
    modd = (x & 1).astype(bool)
    bits = [(x >> b) & 1 for b in range(t+1)]
    f = (p ^ bits[t])[modd]
    s = 1 - 2*f
    L = len(s)
    bias = Fraction(int(s.sum()), L)
    best = None
    for a in range(1, t):
        for b2 in range(a+1, t+1):
            ab = ((bits[a] & bits[b2]) == 1)[modd]
            c1 = int(ab.sum())
            if c1 == 0 or c1 == L: continue
            b1v = Fraction(int(s[ab].sum()), c1)
            b0v = Fraction(int(s[~ab].sum()), L - c1)
            if best is None or abs(b1v) > abs(best[2]):
                best = (a, b2, b1v, b0v)
    a, b2, b1v, b0v = best
    print(f"{t:>3} {str(bias):>12} {('(%d,%d)' % (a, b2)):>11} {str(b1v):>16} {str(b0v):>16}", flush=True)
