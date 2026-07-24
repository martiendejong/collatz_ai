"""
184_lean_export_k13.py
======================
Extend the Lean verification to k=13: gamma = 8619/10000 = 0.8619 (certificate was
certified at 0.86196, so ample margin at the lower rational target).
Same v2 encoding as 183b (flat strings + iterative parser), emitted as a
SECOND lake library `CollatzCert13` in the same project.
Coefficient bounds at gamma = g/d:  A = 2^(-2g/d), B1 = 3^(g/d) 2^(-2g/d),
B3 = 3^(g/d) 2^(-g/d):
   pa^d <= 2^(dQ - 2g)
   p1^d * 2^(2g) <= 3^g * 2^(dQ)
   p3^d * 2^(g)  <= 3^g * 2^(dQ)
"""
import numpy as np
import mpmath as mp
import os

mp.mp.dps = 120
K = 13
N = 3 ** (K - 1)
Nl = N // 3
G, D = 8619, 10000            # gamma = 0.8619 (certified target 0.86196)
Q = 48
S = 34

lam = mp.mpf(2) ** (mp.mpf(G) / D)
alpha = mp.log(3) / mp.log(2)
pa = int(mp.floor(lam ** -2 * 2 ** Q))
p1 = int(mp.floor(lam ** (alpha - 2) * 2 ** Q))
p3 = int(mp.floor(lam ** (alpha - 1) * 2 ** Q))
assert pa ** D <= 2 ** (D * Q - 2 * G)
assert p1 ** D * 2 ** (2 * G) <= 3 ** G * 2 ** (D * Q)
assert p3 ** D * 2 ** G <= 3 ** G * 2 ** (D * Q)

here = os.path.dirname(os.path.abspath(__file__))
v = np.load(os.path.join(here, "certificate_k13.npy")).astype(np.float64)
V = [int(x) for x in np.floor(v * (1 << S)).astype(np.int64)]
assert min(V) > 0 and len(V) == N

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
assert bad == 0, bad
print(f"k=13 pre-check OK: {N:,} inequalities at gamma = {G}/{D}")

root = os.path.join(here, "..", "lean", "CollatzCert")
ddir = os.path.join(root, "CollatzCert13")
os.makedirs(ddir, exist_ok=True)

CHUNK = 20000
nchunks = (N + CHUNK - 1) // CHUNK
for c in range(nchunks):
    lo, hi = c * CHUNK, min((c + 1) * CHUNK, N)
    payload = ",".join(str(V[j]) for j in range(lo, hi)) + ","
    with open(os.path.join(ddir, f"Data{c}.lean"), "w", newline="\n") as f:
        f.write(f'namespace CollatzCert13\n\ndef s{c} : String := "{payload}"\n\n'
                f'end CollatzCert13\n')

imports = "\n".join(f"import CollatzCert13.Data{c}" for c in range(nchunks))
concat = " ++ ".join(f"s{c}" for c in range(nchunks))
main = f"""{imports}

/-!
# Lean-verified feasibility of L_13^NT(2^(8619/10000)) -- pi_1(x) > x^0.8619
Same structure as CollatzCert (k=12); see there for documentation.
-/

namespace CollatzCert13

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

theorem coefA : pa ^ 10000 <= 2 ^ (10000*Q - 17238) := by native_decide
theorem coefB1 : p1 ^ 10000 * 2 ^ 17238 <= 3 ^ 8619 * 2 ^ (10000*Q) := by native_decide
theorem coefB3 : p3 ^ 10000 * 2 ^ 8619 <= 3 ^ 8619 * 2 ^ (10000*Q) := by native_decide
theorem vector_ok : (V.size == N && V.all (fun x => 0 < x)) = true := by native_decide
theorem certificate_feasible : checkAll = true := by native_decide

end CollatzCert13
"""
with open(os.path.join(root, "CollatzCert13.lean"), "w", newline="\n") as f:
    f.write(main)

# register the second library in the lakefile
lakefile = os.path.join(root, "lakefile.toml")
txt = open(lakefile).read()
if "CollatzCert13" not in txt:
    with open(lakefile, "a", newline="\n") as f:
        f.write('\n[[lean_lib]]\nname = "CollatzCert13"\n')
print(f"k=13 Lean library written ({nchunks} chunks)")
