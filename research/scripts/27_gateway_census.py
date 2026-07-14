"""R252-253: the GATEWAY CENSUS.
Every convergent orbit exits through an alternator gate g_y=(4^y-1)/3 (Gateway Law, R250).
Structure: g_y mod 3 cycles 1,2,0 in y -> every third gate (y=3,6,9..: 21, 1365, 87381...)
is a multiple of 3 = a SPRING = no odd predecessors: a door with no hinterland.
Predictions: (1) P(gate y) ~ 4^{-y} for y not=0 mod 3 (subtree mass ~ 1/g);
(2) P(gate y)=0 for y=0 mod 3 except orbits STARTING on the gate's own slide;
(3) w at the door (last macro-step drop) ~ geometric(1/2)."""
import sys, random
from collections import Counter
sys.stdout.reconfigure(encoding="utf-8")
random.seed(252)

gates = {(4**y-1)//3: y for y in range(1, 40)}
N = 200_000
gate_ct, w_door = Counter(), Counter()
for _ in range(N):
    n = random.randrange(3, 1 << 42) | 1
    last, prev = n, None
    m = n
    while m != 1:
        nm = 3*m+1 if m % 2 else m//2
        if nm % 2 and nm != 1: prev, last = last, nm
        m = nm
    y = gates.get(last, -1)
    gate_ct[y] += 1
    if prev is not None and y == gates.get(last):
        # w at the door: halvings between prev (odd) and gate: 3*prev+1 = gate*2^w
        v = 3*prev + 1; w = 0
        while v % 2 == 0: v //= 2; w += 1
        if v == last: w_door[w] += 1

print(f"gateway census, {N:,} random starts in [3, 2^42):")
tot = sum(gate_ct.values())
print(f"{'y':>3} {'gate':>12} {'count':>8} {'P':>10} {'P*4^y':>8}  note")
for y in sorted(k for k in gate_ct if k > 0):
    g = (4**y-1)//3
    p = gate_ct[y]/tot
    note = "SPRING (mult of 3)" if g % 3 == 0 else ""
    print(f"{y:>3} {g:>12,} {gate_ct[y]:>8,} {p:>10.6f} {p*4**y:>8.3f}  {note}")
if gate_ct[-1]: print(f"  NON-GATE last odd: {gate_ct[-1]} <-- would falsify Gateway Law")

print("\nw at the door (drop from penultimate odd into the gate):")
tw = sum(w_door.values())
for w in sorted(w_door)[:8]:
    print(f"  w={w}: {w_door[w]/tw:.4f}   (geometric(1/2) predicts {0.5**w:.4f})")
