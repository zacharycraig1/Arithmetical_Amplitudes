#!/usr/bin/env sage
"""Quick test to get ideal degree."""
F = GF(101)
R = PolynomialRing(F, 'x,y', order='lex')
x, y = R.gens()

# Test degree computation
I = R.ideal([x**3 - 1, y**2 - 1])
print("Test ideal dimension:", I.dimension())
print("Test ideal degree:", I.degree())
print("Expected: dim=0, degree=6")

# Another method for zero-dim ideals
if I.dimension() == 0:
    # Count variety
    K = GF(101**6, 'a')  # extension that splits both
    R2 = PolynomialRing(K, 'x,y', order='lex')
    I2 = R2.ideal([g.change_ring(R2) for g in I.gens()])
    print("Variety size over extension:", len(I2.variety()))
