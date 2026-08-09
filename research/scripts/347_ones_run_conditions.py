# 347: conditions for the appearance of long ones-runs in the combined run-dynamics
# Combined step (course 2.2 form): odd n = a*2^k - 1 (a odd, k trailing ones)
#   -> n' = (a*3^k - 1)/2, then strip remaining zeros to reach the next odd m1.
# Frame: A = a*3^k is odd; A = m1*2^v + 1 with v = v2(A-1).
# So the binary tail of A reads:  1  0^(v-1)  [trailing bits of m1] ...
# and the NEXT ones-run k' is literally the run of ones in A starting at bit v.
# Hypotheses to verify:
#  (H1) reading frame: k' = run of ones of A at offset v, exactly, every step.
#  (H2) for each (k, v, K) there is exactly ONE residue class of a mod 2^(v+K+1)
#       that yields a next run of exactly K ones; hence P(k'=K) ~ 2^-(K+1) style decay.
#  (H3) empirical distribution of k' is geometric ~2^-K and independent of current k.
#  (H4) growth condition: log2-factor of a macro-step is k*log2(3) - (k+v),
#       positive (number grows) iff k >= 2 roughly; tabulate by (k,v).
import random
from math import log2
from collections import defaultdict

def v2(x):
    c = 0
    while x % 2 == 0:
        x //= 2
        c += 1
    return c

def ones_run(x):  # trailing ones of odd x
    c = 0
    while x % 2 == 1:
        c += 1
        x //= 2
    return c

def run_at(A, off):  # run of ones in A starting at bit off
    c = 0
    while (A >> (off + c)) & 1:
        c += 1
    return c

random.seed(347)
freq = defaultdict(int)          # k' frequency
cond = defaultdict(lambda: defaultdict(int))  # k -> k' frequency
drift = defaultdict(list)        # (k,v) -> log2 factors
h1_checked = h2_checked = 0
steps_total = 0

for trial in range(4000):
    n = random.randrange(3, 2**40, 2)
    for _ in range(300):
        if n == 1:
            break
        k = ones_run(n)
        a = (n + 1) >> k
        A = a * 3**k
        v = v2(A - 1)
        m1 = (A - 1) >> v
        # H1: next run read directly from A's bits at offset v
        kp = ones_run(m1) if m1 > 1 else None
        if m1 > 1:
            assert run_at(A, v) == kp, (n, A, v, kp)
            h1_checked += 1
            freq[kp] += 1
            cond[min(k, 6)][kp] += 1
        drift[(min(k, 4), min(v, 4))].append(k * log2(3) - (k + v))
        steps_total += 1
        n = m1

# H2: unique residue class check — for sampled (k, v, K) find all residues of
# a mod 2^(v+K+1) (a odd) whose A = a*3^k has tail 1 0^(v-1) 1^K 0 at the right offsets
import itertools
h2_report = []
for (k, v, K) in [(1, 1, 3), (2, 1, 4), (2, 2, 3), (3, 1, 5), (1, 3, 2), (4, 2, 4)]:
    M = 2 ** (v + K + 1)
    hits = []
    for a in range(1, 2 * M, 2):
        A = a * 3**k
        if v2(A - 1) == v and run_at(A, v) == K:
            hits.append(a % M)
    uniq = sorted(set(hits))
    h2_report.append((k, v, K, len(uniq), uniq[:3]))
    h2_checked += 1

print(f"H1 reading-frame verified on {h1_checked} conversion steps (0 failures)")
print()
print("H2: residues of a mod 2^(v+K+1) giving next-run exactly K (expect exactly 1):")
for k, v, K, cnt, ex in h2_report:
    print(f"  k={k} v={v} K={K}: {cnt} residue class(es) mod 2^{v+K+1}, e.g. {ex}")
print()
tot = sum(freq.values())
print("H3: distribution of next ones-run k' (expect ~2^-(K+1) decay, i.e. half per extra one):")
for K in sorted(freq)[:9]:
    print(f"  k'={K}: {freq[K]/tot:.4f}   (2^-{K+1} = {2**-(K+1):.4f})")
print()
print("H3b: P(k'=K | current k) — rows should match (independence):")
for k in sorted(cond)[:5]:
    row = cond[k]
    t = sum(row.values())
    print(f"  k={k}: " + "  ".join(f"{row.get(K,0)/t:.3f}" for K in range(1, 6)) + f"   (n={t})")
print()
print("H4: mean log2 growth factor per macro-step by (k, v) — positive = number grows:")
for (k, v) in sorted(drift):
    d = drift[(k, v)]
    print(f"  k={k} v={v}: {sum(d)/len(d):+.3f}  (n={len(d)})")
print()
al = [x for d in drift.values() for x in d]
print(f"overall mean drift per macro-step: {sum(al)/len(al):+.4f} log2 (theory: 2*log2(3)-4 = {2*log2(3)-4:+.4f})")
