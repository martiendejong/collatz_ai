import CollatzCert14.Data0
import CollatzCert14.Data1
import CollatzCert14.Data2
import CollatzCert14.Data3
import CollatzCert14.Data4
import CollatzCert14.Data5
import CollatzCert14.Data6
import CollatzCert14.Data7
import CollatzCert14.Data8
import CollatzCert14.Data9
import CollatzCert14.Data10
import CollatzCert14.Data11
import CollatzCert14.Data12
import CollatzCert14.Data13
import CollatzCert14.Data14
import CollatzCert14.Data15
import CollatzCert14.Data16
import CollatzCert14.Data17
import CollatzCert14.Data18
import CollatzCert14.Data19
import CollatzCert14.Data20
import CollatzCert14.Data21
import CollatzCert14.Data22
import CollatzCert14.Data23
import CollatzCert14.Data24
import CollatzCert14.Data25
import CollatzCert14.Data26
import CollatzCert14.Data27
import CollatzCert14.Data28
import CollatzCert14.Data29
import CollatzCert14.Data30
import CollatzCert14.Data31
import CollatzCert14.Data32
import CollatzCert14.Data33
import CollatzCert14.Data34
import CollatzCert14.Data35
import CollatzCert14.Data36
import CollatzCert14.Data37
import CollatzCert14.Data38
import CollatzCert14.Data39
import CollatzCert14.Data40
import CollatzCert14.Data41
import CollatzCert14.Data42
import CollatzCert14.Data43
import CollatzCert14.Data44
import CollatzCert14.Data45
import CollatzCert14.Data46
import CollatzCert14.Data47
import CollatzCert14.Data48
import CollatzCert14.Data49
import CollatzCert14.Data50
import CollatzCert14.Data51
import CollatzCert14.Data52
import CollatzCert14.Data53
import CollatzCert14.Data54
import CollatzCert14.Data55
import CollatzCert14.Data56
import CollatzCert14.Data57
import CollatzCert14.Data58
import CollatzCert14.Data59
import CollatzCert14.Data60
import CollatzCert14.Data61
import CollatzCert14.Data62
import CollatzCert14.Data63
import CollatzCert14.Data64
import CollatzCert14.Data65
import CollatzCert14.Data66
import CollatzCert14.Data67
import CollatzCert14.Data68
import CollatzCert14.Data69
import CollatzCert14.Data70
import CollatzCert14.Data71
import CollatzCert14.Data72
import CollatzCert14.Data73
import CollatzCert14.Data74
import CollatzCert14.Data75
import CollatzCert14.Data76
import CollatzCert14.Data77
import CollatzCert14.Data78
import CollatzCert14.Data79

/-!
# Lean-verified feasibility of L_14^NT(2^(8714/10000)) -- pi_1(x) > x^0.8714
Same structure as CollatzCert (k=12); see there for documentation.
-/

namespace CollatzCert14

def N : Nat := 1594323
def Nl : Nat := 531441
def Q : Nat := 48
def pa : Nat := 84101688550684
def p1 : Nat := 219062569176757
def p3 : Nat := 400761243064555

def dataStr : String := s0 ++ s1 ++ s2 ++ s3 ++ s4 ++ s5 ++ s6 ++ s7 ++ s8 ++ s9 ++ s10 ++ s11 ++ s12 ++ s13 ++ s14 ++ s15 ++ s16 ++ s17 ++ s18 ++ s19 ++ s20 ++ s21 ++ s22 ++ s23 ++ s24 ++ s25 ++ s26 ++ s27 ++ s28 ++ s29 ++ s30 ++ s31 ++ s32 ++ s33 ++ s34 ++ s35 ++ s36 ++ s37 ++ s38 ++ s39 ++ s40 ++ s41 ++ s42 ++ s43 ++ s44 ++ s45 ++ s46 ++ s47 ++ s48 ++ s49 ++ s50 ++ s51 ++ s52 ++ s53 ++ s54 ++ s55 ++ s56 ++ s57 ++ s58 ++ s59 ++ s60 ++ s61 ++ s62 ++ s63 ++ s64 ++ s65 ++ s66 ++ s67 ++ s68 ++ s69 ++ s70 ++ s71 ++ s72 ++ s73 ++ s74 ++ s75 ++ s76 ++ s77 ++ s78 ++ s79

def parseNats (s : String) : Array Nat := Id.run do
  let b := s.toUTF8
  let mut arr : Array Nat := Array.mkEmpty 1594323
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

theorem coefA : pa ^ 10000 <= 2 ^ (10000*Q - 17428) := by native_decide
theorem coefB1 : p1 ^ 10000 * 2 ^ 17428 <= 3 ^ 8714 * 2 ^ (10000*Q) := by native_decide
theorem coefB3 : p3 ^ 10000 * 2 ^ 8714 <= 3 ^ 8714 * 2 ^ (10000*Q) := by native_decide
theorem vector_ok : (V.size == N && V.all (fun x => 0 < x)) = true := by native_decide
theorem certificate_feasible : checkAll = true := by native_decide

end CollatzCert14
