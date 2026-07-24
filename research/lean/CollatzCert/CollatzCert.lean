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

def N : Nat := 177147
def Nl : Nat := 59049
def Q : Nat := 48
def pa : Nat := 86394218259208
def p1 : Nat := 220288582457858
def p3 : Nat := 397621211839769

def dataStr : String := s0 ++ s1 ++ s2 ++ s3 ++ s4 ++ s5 ++ s6 ++ s7 ++ s8

def parseNats (s : String) : Array Nat := Id.run do
  let b := s.toUTF8
  let mut arr : Array Nat := Array.mkEmpty 177147
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
