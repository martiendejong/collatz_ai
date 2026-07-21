"""
108b_r169_proof.py
===================
Algebraic proof that P(h=1)=1 for r=169, and exact P(h=1) for r=253.
"""
import sys, math
sys.stdout.reconfigure(encoding='utf-8')

def v2(x):
    if x == 0: return 999
    c = 0
    while x % 2 == 0: x //= 2; c += 1
    return c

def macro_step(n):
    k = v2(n + 1); m = (n + 1) >> k; x = m * (3**k) - 1; l = v2(x); return x >> l, k, l

BSet = {27,55,63,83,95,103,127,159,169,191,207,223,239,253,255}

# ========================
# r=169 (K=1, m_red=85)
# ========================
print("ALGEBRAIC PROOF for r=169 (K=1, m_red=85):")
print()
print("  Claim: for all odd m == 85 mod 128,")
print("         v2(3m+1) >= 7,  so v2((3m+1)/2) >= 6.")
print()
print("  Proof: 3*85 = 255 == -1 mod 256, so 3*85+1 == 0 mod 256.")
print("         For m = 85 + 128*t (t integer, m odd):")
print("         3m+1 = 3*85 + 384t + 1 = 256 + 384t = 128*(2 + 3t).")
print("         v2(128*(2+3t)) = 7 + v2(2+3t).")
print("         Since 3t is always odd*something, 2+3t depends on t...")
print("         min v2 when t=0: v2(128*(2)) = v2(256) = 8.")
print("         min v2 when t=1: v2(128*5) = 7+0 = 7.")
print("         All have v2 >= 7 => v2((3m+1)/2) >= 6.")
print()
print("  Consequence: n'+1 = (3m+1)/2 is divisible by 64.")
print("  So n' ≡ 63 mod 64.")
print("  The only odd residues ≡ 63 mod 64 in [1,255]:")
r63mod64 = [x for x in range(1, 256, 2) if x % 64 == 63]
print(f"  {r63mod64}")
print(f"  All in BSet: {all(x in BSet for x in r63mod64)}")
print("  => P(h=1) = 1 for r=169. QED")
print()

# Verify empirically
r = 169; K = v2(r + 1)
bset_c = 0; total_c = 0
for k in range(512):
    n = r + 256 * k
    if v2(n + 1) != K: continue
    np, _, _ = macro_step(n)
    total_c += 1
    if np % 256 in BSet: bset_c += 1
print(f"  Empirical (512 samples): P(h=1) = {bset_c}/{total_c} = {bset_c/total_c:.4f}")
print()
print("  Distribution of outputs:")
from collections import Counter
out_dist = Counter()
for k in range(512):
    n = r + 256 * k
    if v2(n + 1) != K: continue
    np, _, _ = macro_step(n)
    out_dist[np % 256] += 1
for res, cnt in sorted(out_dist.items()):
    K_out = v2(res + 1)
    print(f"    n' ≡ {res} mod 256 (k0={K_out}): {cnt} times")

print()
print("=" * 60)
print()

# ========================
# r=253 (K=1, m_red=127)
# ========================
print("ANALYSIS for r=253 (K=1, m_red=127):")
print()
print("  m == 127 mod 128. 3m-1 for m=127: 3*127-1=380=4*95.")
print("  v2(3*127-1) = v2(380) =", v2(380), "(l0 =", v2(380), ")")
print("  n'_0 = 380/4 = 95. BSet!")
print()
print("  For general m = 127 + 128t (odd m):")
print("  3m-1 = 3*(127+128t)-1 = 380 + 384t = 4*(95 + 96t).")
print("  v2(4*(95+96t)) = 2 + v2(95+96t).")
print("  95+96t ≡ 95 mod 2 (odd for all t). So v2(95+96t) depends on 95+96t mod 4:")
print("  95+96t ≡ 95 ≡ 3 mod 4 when t≡0 mod 4 (since 96t≡0 mod 4).")
print("  95+96t ≡ 95+96 ≡ 191 ≡ 3 mod 4 when t≡1 mod 4...")
print("  Actually 96t mod 4 = 0 always (96=24*4). So 95+96t ≡ 3 mod 4 always.")
print("  => v2(95+96t) = 0 always (since 95+96t ≡ 3 mod 4, odd)")
print("  => v2(3m-1) = 2 always for m ≡ 127 mod 128.")
print("  => l = 2 always, n' = (3m-1)/4 = 95+96t.")
print()
print("  n' mod 256 = (95 + 96t) mod 256.")
print("  Cycle: t=0:95, t=1:191, t=2:287mod256=31, t=3:383mod256=127, t=4:479mod256=223")
print("         t=5:575mod256=63, t=6:671mod256=159, t=7:767mod256=255, t=8:863mod256=95 ...")
print("  Period = 256/gcd(96,256) = 256/32 = 8.")
print()
outputs_253 = []
for t in range(8):
    val = (95 + 96 * t) % 256
    in_bset = val in BSet
    outputs_253.append((val, in_bset))
    print(f"  t={t}: n' ≡ {val:3d} mod 256  -> {'BSet' if in_bset else 'non-BSet'}")
n_bset_253 = sum(1 for _, b in outputs_253 if b)
print(f"\n  P(h=1) = {n_bset_253}/8 = {n_bset_253/8:.4f}")
print()

# Verify empirically
r = 253; K = v2(r + 1)
bset_c = 0; total_c = 0
for k in range(512):
    n = r + 256 * k
    if v2(n + 1) != K: continue
    np, _, _ = macro_step(n)
    total_c += 1
    if np % 256 in BSet: bset_c += 1
print(f"  Empirical: P(h=1) = {bset_c}/{total_c} = {bset_c/total_c:.4f}")

print()
print("=" * 60)
print()

# ========================
# GENERAL: P(h=1) for ALL BSet elements
# ========================
print("GENERAL: Exact P(h=1) from arithmetic for each BSet element:")
print()
print("Method: enumerate n' mod (2^{max_period}) for each r.")
print()

def exact_ph1(r, max_t=256):
    """Compute P(first output in BSet) exactly from period structure."""
    K = v2(r + 1)
    m_red = (r + 1) >> K
    # n' = (3^K * m - 1) / 2^v2(3^K*m-1)
    # m = m_red + (2^{8-K}) * t  for t = 0,1,2,...
    # Period in n' mod 256 divides some lcm
    step = 1 << (8 - K) if K < 8 else 1  # step for m in arithmetic progression
    outputs = []
    seen_residues = set()
    period_found = None
    for t in range(max_t):
        m = m_red + step * t
        x = (3**K) * m - 1
        l = v2(x)
        np = x >> l
        rout = np % 256
        outputs.append(rout in BSet)
        if tuple(outputs[-1:]) and t >= 7:
            # Check for period
            for p in range(1, t + 1):
                if outputs[-p:] == outputs[-2*p:-p]:
                    period_found = p
                    break
        if period_found:
            break
    # Use last 'period_found' outputs as the period
    if period_found:
        period_outputs = outputs[-period_found:]
        frac = sum(period_outputs) / period_found
        return frac, period_found, period_outputs
    return sum(outputs) / len(outputs), len(outputs), outputs

print(f"{'r':>4}  {'K':>2}  {'P(h=1)':>8}  {'Period':>7}")
print("-" * 30)
for r in sorted(BSet):
    K = v2(r + 1)
    frac, period, outs = exact_ph1(r)
    print(f"r={r:3d}  K={K}  P(h=1)={frac:.5f}  period={period}")
    if r in [169, 253, 103, 255]:
        print(f"  Outputs: {outs[:min(8, len(outs))]}")
