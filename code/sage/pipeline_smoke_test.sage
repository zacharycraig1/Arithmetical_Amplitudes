"""
Pipeline smoke test using a simple polynomial system.

This tests the eliminant extraction → factorization → cycle type pipeline
without the expensive CHY Groebner computation.

Uses a toy system: x^4 + ax^2 + bx + c = 0 (degree 4, not 24)
"""
from sage.all import *
import json
import os

print("="*60)
print("PIPELINE SMOKE TEST (toy polynomial, not CHY)")
print("="*60)

# Create output directory
os.makedirs("../../data/smoke_test", exist_ok=True)

# Test polynomial: x^4 - 5x^2 + 6 = (x^2-2)(x^2-3)
# Galois group is V4 (Klein four-group)

R_Q = PolynomialRing(QQ, 'x')
x = R_Q.gen()
f_Q = x**4 - 5*x**2 + 6

print(f"\nTest polynomial over Q: {f_Q}")
print(f"Degree: {f_Q.degree()}")
print(f"Discriminant: {f_Q.discriminant()}")

# Factor mod various primes
print("\n" + "-"*60)
print("FACTORIZATION MOD p (Frobenius cycle types)")
print("-"*60)

records = []
for p in [5, 7, 11, 13, 17, 19, 23, 29, 31]:
    F = GF(p)
    R_p = PolynomialRing(F, 'x')
    f_p = R_p(f_Q)
    
    # Check squarefree
    if gcd(f_p, f_p.derivative()) != 1:
        print(f"p={p}: NOT SQUAREFREE (bad reduction)")
        continue
    
    facs = f_p.factor()
    degs = sorted([g.degree() for g, e in facs for _ in range(e)], reverse=True)
    
    rec = {
        "p": int(p),
        "p_mod4": int(p % 4),
        "factor_degrees": [int(d) for d in degs],
        "factorization": str(facs),
    }
    records.append(rec)
    
    print(f"p={p} (mod 4 = {p%4}): degrees={degs}  factors={facs}")

# Save results
results = {
    "polynomial": str(f_Q),
    "degree": int(f_Q.degree()),
    "records": records,
}

outfile = "../../data/smoke_test/toy_frobenius.json"
with open(outfile, "w") as f:
    json.dump(results, f, indent=2)
print(f"\nWrote {outfile}")

# Galois group computation (should be fast for degree 4)
print("\n" + "-"*60)
print("GALOIS GROUP COMPUTATION")
print("-"*60)

try:
    G = f_Q.galois_group()
    print(f"Galois group: {G}")
    print(f"Order: {G.order()}")
    try:
        print(f"Structure: {G.structure_description()}")
    except:
        pass
except Exception as e:
    print(f"Galois group computation failed: {e}")

print("\n" + "="*60)
print("SMOKE TEST COMPLETE")
print("="*60)
print("""
This demonstrates the pipeline:
1. Polynomial definition
2. Factor mod p for multiple primes
3. Extract cycle types from factor degrees
4. Compute Galois group

For CHY, steps 1-3 require Groebner elimination (expensive).
Step 4 requires CRT reconstruction from many primes.
""")
