"""Time the actual CHY eliminant extraction."""
from sage.all import *
import time

# Load kinematics (the load() function handles .sage syntax properly)
load("kinematics.sage")

print("Generating kinematics for seed=42...")
t0 = time.time()
momenta = make_4d_massless_point(7, 42)
mand = momenta_to_mandelstams(momenta)
print(f"Kinematics: {time.time() - t0:.2f}s")

p = 7
F = GF(p)
gauge = (0, 1, -1)

print(f"\nBuilding CHY ideal mod p={p}...")
t1 = time.time()

R = PolynomialRing(F, names=["s7", "s6", "s5", "s4"], order="lex")
s7, s6, s5, s4 = R.gens()

s1, s2, s3 = [F(g) for g in gauge]
sig = {1: s1, 2: s2, 3: s3, 4: s4, 5: s5, 6: s6, 7: s7}

def s_ab(a, b):
    i, j = min(a, b), max(a, b)
    return F(mand[(i, j)])

polys = []
for a in [4, 5, 6, 7]:
    denom_a = R(1)
    for b in [1, 2, 3, 4, 5, 6, 7]:
        if b == a:
            continue
        denom_a *= (sig[a] - sig[b])
    
    num = R(0)
    for b in [1, 2, 3, 4, 5, 6, 7]:
        if b == a:
            continue
        term = s_ab(a, b)
        num += term * (denom_a // (sig[a] - sig[b]))
    
    polys.append(num)

I = Ideal(polys)
print(f"Ideal built: {time.time() - t1:.2f}s")
print(f"Polynomial degrees: {[p.degree() for p in polys]}")

print("\nRunning elimination (this is the slow part)...")
t2 = time.time()
J = I.elimination_ideal([s7, s6, s5])
print(f"Elimination: {time.time() - t2:.2f}s")

gens = J.gens()
univs = [g for g in gens if set(g.variables()) <= set([s4]) and g != 0]
if univs:
    f = min(univs, key=lambda h: h.degree())
    print(f"Eliminant degree: {f.degree()}")
else:
    print("No univariate eliminant found!")
