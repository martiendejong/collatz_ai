import CollatzCert.Data0
import CollatzCert.Data1
import CollatzCert.Data2
import CollatzCert.Data3
import CollatzCert.Data4
import CollatzCert.Data5
import CollatzCert.Data6
import CollatzCert.Data7
import CollatzCert.Data8

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

def N : Nat := 177147
def Nl : Nat := 59049
def Q : Nat := 48
def pa : Nat := 86394218259208
def p1 : Nat := 220288582457858
def p3 : Nat := 397621211839769

def V : Array Nat := data0 ++ data1 ++ data2 ++ data3 ++ data4 ++ data5 ++ data6 ++ data7 ++ data8

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

/-- The certificate satisfies all 177147 inequalities of L_12(2^(213/250)). -/
theorem certificate_feasible : checkAll = true := by native_decide

end CollatzCert
