"""Quick elimination timing test."""
from sage.all import *
import time

p = 7
F = GF(p)
R = PolynomialRing(F, names=['s7','s6','s5','s4'], order='lex')
s7, s6, s5, s4 = R.gens()

print('Building test system...')
t0 = time.time()

# Simple test polynomials (not actual CHY, just for timing)
f1 = s7**3 + s6*s5 + s4
f2 = s6**3 + s7*s4 + s5
f3 = s5**3 + s7*s6 + s4
f4 = s4**3 + s5*s6 + s7

I = Ideal([f1, f2, f3, f4])
print(f'Ideal built: {time.time() - t0:.2f}s')

t1 = time.time()
J = I.elimination_ideal([s7, s6, s5])
print(f'Elimination: {time.time() - t1:.2f}s')
print('Gens:', J.gens())
