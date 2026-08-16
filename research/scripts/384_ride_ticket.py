# 384: RIDE THE 186-BILLION TICKET — construct the unique smallest positive
# integer following the balanced pattern of the first surviving cycle ticket
# (K = 72,057,431,991 / N = 114,208,327,604) for t Syracuse steps.
# Method: required v_i = balanced word; convert to Terras parity string
# ('1' + '0'*(v-1) per step); bit-lift n mod 2^L (bijection).
K = 72057431991
N = 114208327604
p, q = N - K, K
t = 120
vs = [1 + ((i*p)//q - ((i-1)*p)//q) for i in range(1, t+1)]
S = sum(vs)
parity = ''.join('1' + '0'*(v-1) for v in vs) + '1'
L = len(parity)
assert L == S + 1

def Tstep(n):
    return (3*n + 1)//2 if n % 2 else n//2

# bit-lifting: choose bit j so that parity of T^j(r) matches parity[j]
r = 0
for j in range(L):
    n = r
    for _ in range(j):
        n = Tstep(n)
    if n % 2 != int(parity[j]):
        r += 1 << j
    # ensure bit j leaves earlier parities intact (it does: T^i parity depends on r mod 2^(i+1))
n0 = r
print(f"ticket-patroon: K={K:,} N={N:,}; prefix t={t} stappen, {S} gepinde bits")
print(f"kleinste meerijder n0 ({n0.bit_length()} bits):\n  {n0}")

# verify: Syracuse v-sequence matches
from math import log2
n = n0
peaks = [n0]
ok = 0
for i, v_req in enumerate(vs):
    m = 3*n + 1
    v = (m & -m).bit_length() - 1
    if v != v_req:
        break
    n = m >> v
    peaks.append(n)
    ok += 1
print(f"\nverificatie: volgt het ticketpatroon exact {ok}/{t} stappen")
drift = [log2(x/n0) for x in peaks]
print(f"vlakke wiebel (het cyclus-kenmerk): log2(n_i/n0) in [{min(drift):+.3f}, {max(drift):+.3f}]")
print(f"eindpunt van de rit: log2(n_t/n0) = {drift[-1]:+.4f}  (een cyclus zou exact 0.0000 zijn)")

# after the window: the fall
n_end = peaks[-1]
n = n_end
steps_after = 0
while n >= n0 and steps_after < 100000:
    m = 3*n + 1
    v = (m & -m).bit_length() - 1
    n = m >> v
    steps_after += 1
print(f"\nna het venster: zakt onder n0 na {steps_after} extra Syracuse-stappen "
      f"({'gewone sterveling' if steps_after < 100000 else 'nog niet gezakt'})")
print(f"\nvolledige rit zou vergen: {N:,} gepinde bits (~{N/8/2**30:.0f} GB aan getal); "
      f"en het gebalanceerde patroon sluit bewezen NIET (Obs 545)")
