"""
183_lean_export.py
==================
Export the k=12 certificate as an EXACT integer-arithmetic Lean 4 project:
"the first Lean-verified improvement of the Krasikov-Lagarias bound".

Reduction to integers (Obs 349):
  gamma = 213/250 (= 0.852 > 0.84 = published record), lambda = 2^gamma.
  True coefficients: A = lam^-2, B1 = lam^(alpha-2), B3 = lam^(alpha-1),
  alpha = log2(3). Choose dyadic lower bounds a = pa/2^Q <= A etc. Their
  correctness is a pure integer power inequality:
     pa^250 * 2^426 <= 2^(250*Q)            (a^250 <= 2^-426 * 2^...)
     p1^250 * 2^426 <= 3^213 * 2^(250*Q)
     p3^250 * 2^213 <= 3^213 * 2^(250*Q)
  Vector: V[i] = floor(v[i] * 2^S) integers. Feasibility v <= F(v) is implied
  by the integer inequalities (per class i, N = 3^11 = 177147, Nl = N/3):
     D1 (i%3==0): V[i]*2^(2Q) <= pa*2^Q*V[t4 i] + p1*2^Q*Vbar[r1]
     D2 (i%3==1): V[i]*2^Q    <= pa*V[t4 i]
     D3 (i%3==2): V[i]*2^(2Q) <= pa*2^Q*V[t4 i] + p3*2^Q*Vbar[r3]
  ... wait: keep uniform scale 2^(2Q) with pa*2^Q etc. (see check below).
  (Lowering both coefficients and V only ever weakens the RHS relative to the
  real inequality except V on the LHS which is also lowered -- so we simply
  CHECK the integer inequalities exactly; if they pass, the real certificate
  at these dyadic values is feasible, and monotonicity in the coefficients
  gives feasibility at the true (larger) coefficients.)

Pipeline: exact big-int pre-check in Python, then emit a mathlib-free Lean 4
project (core Nat arithmetic + native_decide).
"""
import numpy as np
from math import log2
import mpmath as mp
import os

mp.mp.dps = 120
K = 12
N = 3 ** (K - 1)
Nl = N // 3
GNUM, GDEN = 213, 250            # gamma = 213/250 = 0.852
Q = 48                            # dyadic precision of coefficients
S = 34                            # vector scale 2^S

lam = mp.mpf(2) ** (mp.mpf(GNUM) / GDEN)
alpha = mp.log(3) / mp.log(2)
A_true = lam ** -2
B1_true = lam ** (alpha - 2)
B3_true = lam ** (alpha - 1)

def dyadic_lower(x, q):
    return int(mp.floor(x * (mp.mpf(2) ** q)))

pa = dyadic_lower(A_true, Q)
p1 = dyadic_lower(B1_true, Q)
p3 = dyadic_lower(B3_true, Q)

# exact integer verification of the coefficient bounds
# a = pa/2^Q <= 2^(-2*213/250)  <=>  pa^250 <= 2^(250Q - 426)
c1 = pa ** GDEN <= 2 ** (GDEN * Q - 2 * GNUM)
# b1 = p1/2^Q <= 3^(213/250) * 2^(-426/250) <=> p1^250 * 2^426 <= 3^213 * 2^(250Q)
c2 = p1 ** GDEN * 2 ** (2 * GNUM) <= 3 ** GNUM * 2 ** (GDEN * Q)
# b3 = p3/2^Q <= 3^(213/250) * 2^(-213/250) <=> p3^250 * 2^213 <= 3^213 * 2^(250Q)
c3 = p3 ** GDEN * 2 ** GNUM <= 3 ** GNUM * 2 ** (GDEN * Q)
print(f"coefficient bounds exact-verified: {c1}, {c2}, {c3}")
assert c1 and c2 and c3

v = np.load(os.path.join(os.path.dirname(__file__), "certificate_k12.npy")).astype(np.float64)
assert v.shape[0] == N and v.min() > 0
V = [int(x) for x in np.floor(v * (1 << S)).astype(np.int64)]
assert min(V) > 0
print(f"vector: N={N}, min V={min(V)}, max V={max(V)}")

# exact feasibility check with big ints
i_arr = np.arange(N, dtype=np.int64)
T4 = ((4 * i_arr + 2) % N).astype(np.int64)
def vbar(r):
    return min(V[r], V[r + Nl], V[r + 2 * Nl])

bad = 0
worst = 10.0
for i in range(N):
    t = int(T4[i])
    br = i % 3
    if br == 1:
        lhs = V[i] << Q
        rhs = pa * V[t]
    else:
        lhs = V[i] << (2 * Q)
        if br == 0:
            s = i // 3
            rhs = (pa << Q) * V[t] + (p1 << Q) * vbar((4 * s) % Nl)
        else:
            s = i // 3
            rhs = (pa << Q) * V[t] + (p3 << Q) * vbar((2 * s + 1) % Nl)
    if lhs > rhs:
        bad += 1
    else:
        worst = min(worst, rhs / lhs if lhs else 10.0)
print(f"exact integer feasibility: {N - bad}/{N} pass, min RHS/LHS = {worst:.6f}")
assert bad == 0, f"{bad} failures -- increase S or lower gamma"

# ---------------- emit Lean project ----------------------------------------
root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "lean", "CollatzCert"))
os.makedirs(os.path.join(root, "CollatzCert"), exist_ok=True)

with open(os.path.join(root, "lean-toolchain"), "w", newline="\n") as f:
    f.write("leanprover/lean4:v4.15.0\n")

with open(os.path.join(root, "lakefile.toml"), "w", newline="\n") as f:
    f.write('name = "CollatzCert"\ndefaultTargets = ["CollatzCert"]\n\n'
            '[[lean_lib]]\nname = "CollatzCert"\n')

CHUNK = 20000
nchunks = (N + CHUNK - 1) // CHUNK
for c in range(nchunks):
    lo, hi = c * CHUNK, min((c + 1) * CHUNK, N)
    with open(os.path.join(root, "CollatzCert", f"Data{c}.lean"), "w", newline="\n") as f:
        f.write(f"namespace CollatzCert\n\ndef data{c} : Array Nat := #[\n")
        f.write(",\n".join(
            "  " + ", ".join(str(V[j]) for j in range(a, min(a + 20, hi)))
            for a in range(lo, hi, 20)))
        f.write("\n]\n\nend CollatzCert\n")

imports = "\n".join(f"import CollatzCert.Data{c}" for c in range(nchunks))
appends = " ++ ".join(f"data{c}" for c in range(nchunks))
main = f"""{imports}

/-!
# A Lean-verified feasibility certificate for the Krasikov-Lagarias system L_12
For gamma = 213/250 = 0.852 (> 0.84, the published record exponent of
Krasikov-Lagarias, Acta Arith. 109 (2003) 237-258), lambda = 2^gamma, the
inequalities below are an exact integer-arithmetic witness that the system
L_12^NT(lambda) (Prop. 2.1 / (2.7)-(2.14) of loc. cit.) is feasible.
By Theorem 2.2 of loc. cit. (not formalized here; classical, peer-reviewed),
feasibility implies pi_1(x) > x^0.852 for all large x.

Structure: coefficients A = lambda^-2, B1 = lambda^(log2 3 - 2),
B3 = lambda^(log2 3 - 1) are lower-bounded by dyadics pa/2^Q, p1/2^Q, p3/2^Q
(theorems coefA/coefB1/coefB3: pure integer power inequalities), and the
177147-entry certificate vector V satisfies the per-class inequalities
(theorem certificate_feasible). All checks are by native_decide over Nat.
-/

namespace CollatzCert

def N : Nat := {N}
def Nl : Nat := {Nl}
def Q : Nat := {Q}
def pa : Nat := {pa}
def p1 : Nat := {p1}
def p3 : Nat := {p3}

def V : Array Nat := {appends}

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

/-- pa/2^Q <= 2^(-426/250): pure integer power inequality. -/
theorem coefA : pa ^ 250 <= 2 ^ (250*Q - 426) := by native_decide

/-- p1/2^Q <= 3^(213/250) * 2^(-426/250). -/
theorem coefB1 : p1 ^ 250 * 2 ^ 426 <= 3 ^ 213 * 2 ^ (250*Q) := by native_decide

/-- p3/2^Q <= 3^(213/250) * 2^(-213/250). -/
theorem coefB3 : p3 ^ 250 * 2 ^ 213 <= 3 ^ 213 * 2 ^ (250*Q) := by native_decide

/-- Vector is fully positive and of the right length. -/
theorem vector_ok : V.size = N && V.all (fun x => 0 < x) := by native_decide

/-- The certificate satisfies all {N} inequalities of L_12(2^(213/250)). -/
theorem certificate_feasible : checkAll = true := by native_decide

end CollatzCert
"""
with open(os.path.join(root, "CollatzCert.lean"), "w", newline="\n") as f:
    f.write(main)

sizes = sum(os.path.getsize(os.path.join(root, "CollatzCert", x))
            for x in os.listdir(os.path.join(root, "CollatzCert")))
print(f"Lean project written to {root}: {nchunks} data chunks, {sizes/1e6:.1f} MB data")
