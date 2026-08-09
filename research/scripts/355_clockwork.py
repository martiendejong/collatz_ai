# 355: the clockwork — Martien's reformulation.
# Per odd step the factor is (3 + 1/n)/2 = (3/2)*(1 + 1/(3n)): the clock hand
# advances by log2(3) (irrational) and the +1 supplies a tiny extra push
# log2(1+1/(3n)) ~ 1/(3n ln2). A cycle = the clock returning EXACTLY to zero:
#   N - K*log2(3) = sum_i log2(1 + 1/(3 n_i))          (all pushes same sign)
# After K ticks the clock misses integer-N by:
#   delta+(K) = ceil(K x) - K x   (needs POSITIVE elements to close)
#   delta-(K) = K x - floor(K x)  (needs NEGATIVE elements to close)
# Closing requires an element scale  n~ = K / (3 ln2 delta): computable!
# So the "spectrum" of possible cycles is a table (K, n~) — and the known
# cycles (1, -1, -5, -17) should appear at the first entries.
from decimal import Decimal, getcontext
getcontext().prec = 130
x = Decimal(3).ln() / Decimal(2).ln()
ln2 = Decimal(2).ln()

def spectrum_row(K):
    Kx = K * x
    fl = int(Kx)
    dminus = Kx - fl
    dplus = 1 - dminus
    npos = K / (3 * ln2 * dplus)
    nneg = K / (3 * ln2 * dminus)
    return dplus, npos, dminus, nneg

known = {1: "n=1 (trivial 1-4-2-1) / n=-1", 2: "n=-5 cycle", 7: "n=-17 cycle"}
print("the clockwork spectrum, small K (exhaustive):")
print("  K   miss+      scale n~+     miss-      scale n~-    known cycle at this tick?")
for K in range(1, 31):
    dp, npos, dm, nneg = spectrum_row(K)
    mark = known.get(K, "")
    print(f" {K:>2}  {float(dp):8.5f}  {float(npos):12.2f}  {float(dm):8.5f}  {float(nneg):12.2f}   {mark}")
print()
# exact closure checks for the known cycles
def check(elems, N):
    prod = Decimal(1)
    for n in elems:
        prod *= (1 + Decimal(1) / (3 * n))
    lhs = Decimal(2) ** N
    rhs = Decimal(3) ** len(elems) * prod
    return abs(lhs - rhs) < Decimal("1e-90")
print("exact closure checks:")
print(f"  trivial (n=1, K=1, N=2):    {check([1], 2)}")
print(f"  -1 cycle (K=1, N=1):        {check([-1], 1)}")
print(f"  -5 cycle (K=2, N=3):        {check([-5, -7], 3)}")
print(f"  -17 cycle (K=7, N=11):      {check([-17, -25, -37, -55, -41, -61, -91], 11)}")
print()
# the ladder of near-closures: convergents + semiconvergents of x, positive side,
# with the element scale each would require; first scale beyond 2^68 = first
# candidate not excluded by verification
cf_x = x
a0 = int(cf_x)
frac = cf_x - a0
h0, h1, k0, k1 = 1, a0, 0, 1
convs = [(a0, 1)]
while k1 < 10**13:
    y = 1 / frac
    ai = int(y)
    frac = y - ai
    h0, h1 = h1, ai * h1 + h0
    k0, k1 = k1, ai * k1 + k0
    convs.append((h1, k1))
cands = set()
for (p1, q1), (p2, q2) in zip(convs, convs[1:]):
    for t in range(0, 200):
        q = q1 + t * q2
        if q > 2 * 10**11:
            break
        cands.add(q)
V = Decimal(2) ** 68
print("positive-side near-closures (the clock's best ticks) and required element scale:")
print("  K                 miss+          required n~+        vs verified 2^68")
best_shown = 0
for K in sorted(cands):
    dp, npos, dm, nneg = spectrum_row(K)
    if float(npos) < 3 * best_shown + 3:
        continue
    best_shown = float(npos)
    status = "EXCLUDED by verification" if npos < V else "<-- FIRST SURVIVOR"
    print(f"  {K:>13,}  {float(dp):12.3E}  {float(npos):16.3E}   {status}")
    if npos >= V:
        break
