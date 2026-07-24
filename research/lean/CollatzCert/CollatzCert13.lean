import CollatzCert13.Data0
import CollatzCert13.Data1
import CollatzCert13.Data2
import CollatzCert13.Data3
import CollatzCert13.Data4
import CollatzCert13.Data5
import CollatzCert13.Data6
import CollatzCert13.Data7
import CollatzCert13.Data8
import CollatzCert13.Data9
import CollatzCert13.Data10
import CollatzCert13.Data11
import CollatzCert13.Data12
import CollatzCert13.Data13
import CollatzCert13.Data14
import CollatzCert13.Data15
import CollatzCert13.Data16
import CollatzCert13.Data17
import CollatzCert13.Data18
import CollatzCert13.Data19
import CollatzCert13.Data20
import CollatzCert13.Data21
import CollatzCert13.Data22
import CollatzCert13.Data23
import CollatzCert13.Data24
import CollatzCert13.Data25
import CollatzCert13.Data26

/-!
# Lean-verified feasibility of L_13^NT(2^(8619/10000)) -- pi_1(x) > x^0.8619
Same structure as CollatzCert (k=12); see there for documentation.
-/

namespace CollatzCert13

def N : Nat := 531441
def Nl : Nat := 177147
def Q : Nat := 48
def pa : Nat := 85216616239011
def p1 : Nat := 219662081583242
def p3 : Nat := 399220512696425

def dataStr : String := s0 ++ s1 ++ s2 ++ s3 ++ s4 ++ s5 ++ s6 ++ s7 ++ s8 ++ s9 ++ s10 ++ s11 ++ s12 ++ s13 ++ s14 ++ s15 ++ s16 ++ s17 ++ s18 ++ s19 ++ s20 ++ s21 ++ s22 ++ s23 ++ s24 ++ s25 ++ s26

def parseNats (s : String) : Array Nat := Id.run do
  let b := s.toUTF8
  let mut arr : Array Nat := Array.mkEmpty 531441
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
