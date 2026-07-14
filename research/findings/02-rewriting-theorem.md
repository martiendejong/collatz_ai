# F2 — The Rewriting Theorem: binary a·1^k → ternary (a−1)·2^k

**Status: PROVEN (one-line proof) + verified for all a < 400, k ≤ 11.**

## Statement
The segment map sends the binary word ⟨bits of a⟩⟨k ones⟩ to the ternary word ⟨trits of a−1⟩⟨k twos⟩:

> starter a·2^k − 1  ⟶ (segment) ⟶  endpoint a·3^k − 1 = (a−1)·3^k + (3^k − 1)

so the endpoint's ternary representation is literally the digits of (a−1) followed by k twos. ∎

## Examples (from the original chat mystery — now explained)
| starter | binary | endpoint | ternary |
|---------|--------|----------|---------|
| 7 | 111 | 26 | **222** |
| 1 | 1 | 2 | **2** |
| 9 | 1001 | 14 | 11**2** |
| 39 | 100111 | 134 | 11**222** |
| 55 | 110111 | 188 | 20**222** |

The 111₂ → 222₃ observation that started the whole framework is the case a=1, k=3. The trivial cycle is the fixed point of this rewriting system: word "1" → word "2" → (halve) → word "1".

## Complexity localization
The segment (binary→ternary direction) is **trivially exact** — prefix decrements by 1, run converts 1s→2s. ALL of the difficulty of the Collatz problem is concentrated in the return direction: the **halving cascade** x → x/2^w read in ternary, i.e. division by 2 in base 3, whose borrow chains are equidistribution-hard. The problem is not "what does ×3 do to binary" (answered above) but "what does ÷2 do to ternary" (equivalent to the full conjecture).

Related: [[01-family-pair-coordinates]], [[06-ternary-streaks]]
