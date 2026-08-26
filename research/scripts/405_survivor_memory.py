"""Obs 613: Rule 30 lens, task S4 -- branch-point census on the needle.

Question. Thm 20/22: at density level the machine stores nothing
(refills exactly fair and independent). Rule 30's left half teaches
that seed information can still survive at RARE configurations (their
branch points), found only by conditioning on the exceptional set.
Our exceptional set is the needle: orbits that never drop below their
start. Does the survivor ensemble carry serial memory in its k-stream,
beyond what survival-conditioning itself induces through altitude?

Design (controls doctrine).
  - Real: all odd n in two disjoint seed windows; macro-steps
    k = v2(3n+1); survivor prefix = steps while n_i >= n_0.
  - Null: iid k with P(k) = 2^-k, height walk h += log2(3) - k,
    survival = h >= 0 throughout; same estimator, same sample sizes.
    Conditioning a walk to stay positive induces memory THROUGH the
    height, so the null is the honest baseline, not zero.
  - Statistic: plug-in mutual information I(k_i; k_{i+d}) on survivor
    steps, (a) unconditional, (b) conditioned on altitude bin
    (weighted average over bins). k capped at 6, altitude in 6 bins.
  - Permutation floor: shuffle the second coordinate within each
    stratum; excess = MI - MI_perm is the reported signal.
  - Two windows guard against the orbit-sampling contamination trap
    (merged orbits share segments; R187 artifact).

Verdict rule. If Collatz excess matches the null excess in (a) and
both conditional excesses (b) are at floor, then altitude screens off
all memory: the branch-point census on the needle is EMPTY -- the
machine remembers nothing but its own altitude (pointwise sharpening
of Thm 22). Any stable surplus over the null = a lead.
"""
import math
import random

KCAP = 6      # k values 1..5, >=6 pooled
HBINS = 6     # altitude bins
DEPTH = 60    # max survivor depth followed
LAGS = (1, 2, 4)
LOG23 = math.log2(3)

def v2(x: int) -> int:
    v = 0
    while x % 2 == 0:
        x //= 2
        v += 1
    return v

def survivor_prefix_real(n0: int):
    """Yield (k_i, h_i) along the never-drop prefix, h = log2(n_i/n_0)."""
    n = n0
    out = []
    for _ in range(DEPTH):
        k = v2(3 * n + 1)
        h = math.log2(n / n0)
        out.append((min(k, KCAP), h))
        n = (3 * n + 1) >> k
        if n < n0:
            break
    return out

def survivor_prefix_null(rng):
    h = 0.0
    out = []
    for _ in range(DEPTH):
        k = 1
        while rng.random() < 0.5:
            k += 1
        out.append((min(k, KCAP), h))
        h += LOG23 - k
        if h < 0:
            break
    return out

def collect(prefixes, lag):
    """(k_i, k_{i+lag}, hbin_i) triples over surviving steps."""
    triples = []
    for pref in prefixes:
        for i in range(len(pref) - lag):
            k1, h = pref[i]
            k2, _ = pref[i + lag]
            hb = min(HBINS - 1, max(0, int(h / 2.0)))  # 2-bit altitude bins
            triples.append((k1, k2, hb))
    return triples

def plugin_mi(pairs):
    from collections import Counter
    n = len(pairs)
    cj = Counter(pairs); c1 = Counter(a for a, _ in pairs); c2 = Counter(b for _, b in pairs)
    mi = 0.0
    for (a, b), c in cj.items():
        p = c / n
        mi += p * math.log2(p * n * n / (c1[a] * c2[b]))
    return mi

def mi_with_floor(pairs, rng):
    mi = plugin_mi(pairs)
    shuf = [b for _, b in pairs]
    rng.shuffle(shuf)
    floor = plugin_mi(list(zip([a for a, _ in pairs], shuf)))
    return mi, floor

def cond_mi_with_floor(triples, rng):
    """Altitude-weighted conditional MI and its permutation floor."""
    from collections import defaultdict
    strata = defaultdict(list)
    for k1, k2, hb in triples:
        strata[hb].append((k1, k2))
    n = len(triples)
    mi = fl = 0.0
    for pairs in strata.values():
        if len(pairs) < 50:
            continue
        m, f = mi_with_floor(pairs, rng)
        w = len(pairs) / n
        mi += w * m
        fl += w * f
    return mi, fl

def run_window(label, seeds_real=None, n_null=0, rng=None):
    if seeds_real is not None:
        prefixes = [survivor_prefix_real(n) for n in seeds_real]
    else:
        prefixes = [survivor_prefix_null(rng) for _ in range(n_null)]
    # keep only orbits that survive at least a few steps
    prefixes = [p for p in prefixes if len(p) >= 8]
    print(f"{label}: {len(prefixes)} orbits with survivor depth >= 8")
    for lag in LAGS:
        tr = collect(prefixes, lag)
        mi, fl = mi_with_floor([(a, b) for a, b, _ in tr], rng)
        cmi, cfl = cond_mi_with_floor(tr, rng)
        print(f"  lag {lag}: N={len(tr):7d}  MI={mi:.5f} (floor {fl:.5f}, "
              f"excess {mi - fl:+.5f})  MI|h={cmi:.5f} (floor {cfl:.5f}, "
              f"excess {cmi - cfl:+.5f})")

if __name__ == "__main__":
    rng = random.Random(2026)
    W = 400_000
    for base in (1_000_001, 30_000_001):  # two disjoint seed windows
        seeds = range(base, base + 2 * W, 2)
        run_window(f"REAL window {base}", seeds_real=seeds, rng=rng)
    run_window("NULL fair-coin walk", n_null=W, rng=rng)
    run_window("NULL fair-coin walk (replicate)", n_null=W, rng=rng)
