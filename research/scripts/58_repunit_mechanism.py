"""R531-540: THE BINARY MECHANISM OF 3^k (the repunit's image).
Repunit 2^k-1 (k ones) --ladder burn--> 3^k-1 = '222..2' ternary. In binary?
(1) worked list of 3^k-1 in binary;
(2) TRAILING-BIT CLOCKS: bit i of 3^k is periodic in k (3 is a 2-adic unit,
    ord of 3 mod 2^j = 2^(j-2)) -- measure each bit's period = the odometer law;
(3) LENGTH INCREMENTS: bitlength grows by the Sturmian word of log2(3):
    pattern of 2s and 1s = irrational rotation -- THE non-mechanism upstairs;
(4) ones-fraction -> 1/2 (never all-ones again)."""
import sys, math
sys.stdout.reconfigure(encoding="utf-8")

print("(1) 3^k - 1 in binary (the repunit's landing):")
for k in range(1, 21):
    v = 3**k - 1
    print(f"  k={k:>2}: {bin(v)[2:]:>33}  len={v.bit_length()} ones={bin(v).count('1')}")

print("\n(2) trailing-bit clocks of 3^k (bit i as function of k, k=0..32):")
for i in range(0, 8):
    bits = "".join(str((3**k >> i) & 1) for k in range(33))
    # detect period
    per = None
    for p in (1, 2, 4, 8, 16, 32):
        if all(bits[t] == bits[t % p] for t in range(len(bits))):
            per = p; break
    print(f"  bit {i}: {bits}  period {per}")

print("\n(3) bit-length increments of 3^k (Sturmian word of log2 3):")
incs = "".join(str((3**(k+1)).bit_length() - (3**k).bit_length()) for k in range(40))
print(f"  {incs}")
frac = incs.count('2') / len(incs)
print(f"  fraction of 2s: {frac:.3f}  (theory: log2(3)-1 = {math.log2(3)-1:.3f})")

print("\n(4) ones-fraction of 3^k-1:")
for k in (5, 10, 20, 40, 80, 160):
    v = 3**k - 1
    print(f"  k={k:>3}: ones/len = {bin(v).count('1')}/{v.bit_length()} = {bin(v).count('1')/v.bit_length():.3f}")
