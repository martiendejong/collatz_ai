"""
183b_lean_export_v2.py
======================
V2 of the Lean export: the v1 array-literal encoding stack-overflowed native
evaluation (deeply nested initializers, exit 0xC0000409). V2 stores the vector
as flat STRING literals parsed by an iterative ByteArray loop -- no recursion
anywhere. Same verified constants as 183 (re-derived + re-checked here).
"""
import numpy as np
import mpmath as mp
import os

mp.mp.dps = 120
K = 12
N = 3 ** (K - 1)
Nl = N // 3
GNUM, GDEN = 213, 250
Q = 48
S = 34

lam = mp.mpf(2) ** (mp.mpf(GNUM) / GDEN)
alpha = mp.log(3) / mp.log(2)
pa = int(mp.floor(lam ** -2 * 2 ** Q))
p1 = int(mp.floor(lam ** (alpha - 2) * 2 ** Q))
p3 = int(mp.floor(lam ** (alpha - 1) * 2 ** Q))
assert pa ** GDEN <= 2 ** (GDEN * Q - 2 * GNUM)
assert p1 ** GDEN * 2 ** (2 * GNUM) <= 3 ** GNUM * 2 ** (GDEN * Q)
assert p3 ** GDEN * 2 ** GNUM <= 3 ** GNUM * 2 ** (GDEN * Q)

here = os.path.dirname(os.path.abspath(__file__))
v = np.load(os.path.join(here, "certificate_k12.npy")).astype(np.float64)
V = [int(x) for x in np.floor(v * (1 << S)).astype(np.int64)]
assert min(V) > 0 and len(V) == N

# exact big-int feasibility re-check (same as 183)
bad = 0
for i in range(N):
    t = (4 * i + 2) % N
    br = i % 3
    if br == 1:
        ok = V[i] << Q <= pa * V[t]
    else:
        s = i // 3
        r = (4 * s) % Nl if br == 0 else (2 * s + 1) % Nl
        pb = p1 if br == 0 else p3
        vb = min(V[r], V[r + Nl], V[r + 2 * Nl])
        ok = V[i] << (2 * Q) <= (pa << Q) * V[t] + (pb << Q) * vb
    bad += (not ok)
assert bad == 0
print(f"pre-check OK: {N:,} inequalities, coefficients verified")

root = os.path.join(here, "..", "lean", "CollatzCert")
ddir = os.path.join(root, "CollatzCert")
for f in os.listdir(ddir):
    if f.startswith("Data"):
        os.remove(os.path.join(ddir, f))

CHUNK = 20000
nchunks = (N + CHUNK - 1) // CHUNK
for c in range(nchunks):
    lo, hi = c * CHUNK, min((c + 1) * CHUNK, N)
    payload = ",".join(str(V[j]) for j in range(lo, hi)) + ","
    with open(os.path.join(ddir, f"Data{c}.lean"), "w", newline="\n") as f:
        f.write(f'namespace CollatzCert\n\ndef s{c} : String := "{payload}"\n\nend CollatzCert\n')

imports = "\n".join(f"import CollatzCert.Data{c}" for c in range(nchunks))
concat = " ++ ".join(f"s{c}" for c in range(nchunks))
main = f"""{imports}

/-!
# A Lean-verified feasibility certificate for the Krasikov-Lagarias system L_12
gamma = 213/250 = 0.852 > 0.84 (the published record exponent, Krasikov-
Lagarias, Acta Arith. 109 (2003) 237-258; still the record per arXiv:2512.13760).
Feasibility of L_12^NT(2^gamma) is witnessed in exact integer arithmetic:
coefficient lower bounds via integer power inequalities (coefA/coefB1/coefB3),
and the 177,147 per-class inequalities via `certificate_feasible`.
By Theorem 2.2 of loc. cit. (classical, cited, not formalized) this yields
pi_1(x) > x^0.852 for all large x. Data is stored as flat string literals and
parsed by an iterative byte loop (no deep recursion).
-/

namespace CollatzCert

def N : Nat := {N}
def Nl : Nat := {Nl}
def Q : Nat := {Q}
def pa : Nat := {pa}
def p1 : Nat := {p1}
def p3 : Nat := {p3}

def dataStr : String := {concat}

def parseNats (s : String) : Array Nat := Id.run do
  let b := s.toUTF8
  let mut arr : Array Nat := Array.mkEmpty {N}
  let mut cur : Nat := 0
  let mut has : Bool := false
  for i in [0:b.size] do
    let c := b.get! i |>.toNat
    if 48 <= c && c <= 57 then
      cur := cur * 10 + (c - 48)
      has := true
    else
      if has then arr := arr.push cur
      cur := 0
      has := false
  if has then arr := arr.push cur
  return arr

def V : Array Nat := parseNats dataStr

def vbar (r : Nat) : Nat :=
  min (V.getD r 0) (min (V.getD (r + Nl) 0) (V.getD (r + 2*Nl) 0))

def checkOne (i : Nat) : Bool :=
  let t := (4*i + 2) % N
  match i % 3 with
  | 1 => V.getD i 0 <<< Q <= pa * V.getD t 0
  | 0 => V.getD i 0 <<< (2*Q) <= (pa <<< Q) * V.getD t 0 + (p1 <<< Q) * vbar ((4*(i/3)) % Nl)
  | _ => V.getD i 0 <<< (2*Q) <= (pa <<< Q) * V.getD t 0 + (p3 <<< Q) * vbar ((2*(i/3)+1) % Nl)

def checkAll : Bool := Id.run do
  let mut ok := true
  for i in [0:N] do
    ok := ok && checkOne i
  return ok

/-- pa/2^Q <= lambda^-2 = 2^(-426/250): integer power inequality. -/
theorem coefA : pa ^ 250 <= 2 ^ (250*Q - 426) := by native_decide

/-- p1/2^Q <= lambda^(log2 3 - 2) = 3^(213/250) * 2^(-426/250). -/
theorem coefB1 : p1 ^ 250 * 2 ^ 426 <= 3 ^ 213 * 2 ^ (250*Q) := by native_decide

/-- p3/2^Q <= lambda^(log2 3 - 1) = 3^(213/250) * 2^(-213/250). -/
theorem coefB3 : p3 ^ 250 * 2 ^ 213 <= 3 ^ 213 * 2 ^ (250*Q) := by native_decide

/-- The parsed vector has the right length and is fully positive. -/
theorem vector_ok : (V.size == N && V.all (fun x => 0 < x)) = true := by native_decide

/-- All 177,147 inequalities of L_12^NT(2^(213/250)) hold for the certificate. -/
theorem certificate_feasible : checkAll = true := by native_decide

end CollatzCert
"""
with open(os.path.join(root, "CollatzCert.lean"), "w", newline="\n") as f:
    f.write(main)
print(f"v2 Lean project written ({nchunks} string chunks)")
