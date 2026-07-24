"""
185_k21_memmap.py
=================
k=21 Perron iteration for L_21^NT(lambda), memmapped on disk (HDD-friendly:
ALL access strictly sequential via the affine window structure).

  N = 3^20 = 3,486,784,401 classes; Nl = 3^19.
  Backbone gather v[(4i+2) mod N]: consecutive output chunks read consecutive
  input windows -- the sweep reads v exactly 4x sequentially (wrapping 3x).
  Feeds r1 = 4s mod Nl (stride-4 slice of a window), r3 = 2s+1 (stride-2).

Runbook lessons applied (RUNBOOK_k21.md): checkpoint state every sweep (files
persist; state.json carries parity+norm), chunked everything, warm start by
3x-tiling the k=20 polished vector (period 3^19 = Nl), margin tracked cheaply
every RATIO_EVERY sweeps, target lambda = 1.890 (pre-registered prediction
gamma(21) ~ 0.918-0.919, PREDICTIONS.md #1).

Usage: python 185_k21_memmap.py [sweeps] [lam_num]   (defaults 40, 1890)
State/files in research/k21/ on E: (1.6TB free).
"""
import numpy as np
import json, os, sys, time
from math import log2

K = 21
N = 3 ** (K - 1)          # 3,486,784,401
Nl = N // 3               # 1,162,261,467
CH = 3 ** 13              # 1,594,323 outputs per chunk
NCH = N // CH             # 2187
ALPHA = log2(3.0)

HERE = os.path.dirname(os.path.abspath(__file__))
DIR = os.path.abspath(os.path.join(HERE, "..", "k21"))
os.makedirs(DIR, exist_ok=True)
STATE = os.path.join(DIR, "state.json")

SWEEPS = int(sys.argv[1]) if len(sys.argv) > 1 else 40
LAM_NUM = int(sys.argv[2]) if len(sys.argv) > 2 else 1890
lam = LAM_NUM / 1000.0
A = np.float32(lam ** -2)
B1 = np.float32(lam ** (ALPHA - 2))
B3 = np.float32(lam ** (ALPHA - 1))
RATIO_EVERY = 5

def open_vec(name, create=False):
    p = os.path.join(DIR, name)
    if create and not os.path.exists(p):
        return np.lib.format.open_memmap(p, mode="w+", dtype=np.float32, shape=(N,))
    return np.lib.format.open_memmap(p, mode="r+")

def read_wrap(vec, start, length):
    """Sequential read of [start, start+length) mod N (wraps at most once)."""
    start %= N
    end = start + length
    if end <= N:
        return np.asarray(vec[start:end])
    a = np.asarray(vec[start:N])
    b = np.asarray(vec[:end - N])
    return np.concatenate([a, b])

def read_wrap_l(vec, start, length):
    start %= Nl
    end = start + length
    if end <= Nl:
        return np.asarray(vec[start:end])
    return np.concatenate([np.asarray(vec[start:Nl]), np.asarray(vec[:end - Nl])])

# ---------------- init / resume ---------------------------------------------
if os.path.exists(STATE):
    st = json.load(open(STATE))
    print(f"resuming at sweep {st['sweep']}, norm {st['norm']:.6g}", flush=True)
else:
    print("initializing: warm start = 3x tile of k20 polished vector", flush=True)
    src = np.lib.format.open_memmap(
        os.path.join(HERE, "..", "certificates", "k20_polished.npy"), mode="r")
    assert src.shape[0] == Nl
    va = open_vec("va.npy", create=True)
    for lo in range(0, Nl, CH * 3):
        hi = min(lo + CH * 3, Nl)
        block = np.asarray(src[lo:hi], dtype=np.float32)
        va[lo:hi] = block
        va[Nl + lo:Nl + hi] = block
        va[2 * Nl + lo:2 * Nl + hi] = block
    va.flush(); del va, src
    st = {"sweep": 0, "cur": "va.npy", "oth": "vb.npy", "norm": 1.0}
    open_vec("vb.npy", create=True)
    json.dump(st, open(STATE, "w"))
    print("warm start written", flush=True)

cbar_path = os.path.join(DIR, "cbar.npy")
if os.path.exists(cbar_path):
    cbar = np.lib.format.open_memmap(cbar_path, mode="r+")
else:
    cbar = np.lib.format.open_memmap(cbar_path, mode="w+", dtype=np.float32,
                                     shape=(Nl,))

print(f"k=21: lambda = {lam} (gamma = {log2(lam):.5f}), {SWEEPS} sweeps, "
      f"chunk {CH:,} x {NCH}", flush=True)

for sweep in range(st["sweep"], SWEEPS):
    t0 = time.time()
    v = open_vec(st["cur"])
    w = open_vec(st["oth"])
    inv_norm = np.float32(1.0 / st["norm"])

    # pass 1: cbar[r] = min over the three lifts (normalized)
    for lo in range(0, Nl, CH):
        hi = min(lo + CH, Nl)
        m0 = np.asarray(v[lo:hi])
        m1 = np.asarray(v[Nl + lo:Nl + hi])
        m2 = np.asarray(v[2 * Nl + lo:2 * Nl + hi])
        cbar[lo:hi] = np.minimum(np.minimum(m0, m1), m2) * inv_norm

    # pass 2: w = F(v/norm), tracking max and (periodically) min ratio
    wmax = 0.0
    rmin = np.inf
    track = (sweep % RATIO_EVERY == 0) or (sweep == SWEEPS - 1)
    for c in range(NCH):
        lo = c * CH
        hi = lo + CH
        win = read_wrap(v, 4 * lo + 2, 4 * CH)      # backbone window
        out = win[::4].astype(np.float32) * (A * inv_norm)
        # D1: i = 3j (rel. offsets 0,3,6,...): s = lo/3 + j, r1 = 4s mod Nl
        w1 = read_wrap_l(cbar, 4 * (lo // 3), 4 * (CH // 3))[::4]
        out[0::3] += B1 * w1
        # D3: i = 3j+2: s = (lo+2)/3? s = i//3 = lo//3 + j, r3 = 2s+1 mod Nl
        w3 = read_wrap_l(cbar, 2 * (lo // 3) + 1, 2 * (CH // 3))[::2]
        out[2::3] += B3 * w3
        w[lo:hi] = out
        m = float(out.max())
        if m > wmax:
            wmax = m
        if track:
            vv = np.asarray(v[lo:hi]) * inv_norm
            r = float((out / vv).min())
            if r < rmin:
                rmin = r
    st = {"sweep": sweep + 1, "cur": st["oth"], "oth": st["cur"], "norm": wmax}
    json.dump(st, open(STATE, "w"))
    msg = f"sweep {sweep + 1}/{SWEEPS}  growth={wmax:.6f}"
    if track:
        msg += f"  min F(v)/v = {rmin:.6f}"
    msg += f"  [{time.time() - t0:.0f}s]"
    print(msg, flush=True)

print("done. growth -> rho(lambda); feasible iff rho >= 1 "
      "(then floor + exact verify, script 186).", flush=True)
