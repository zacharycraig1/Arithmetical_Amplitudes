#!/usr/bin/env sage
"""
Simple generic kinematics using known parameterization
"""

import sys
import time

print("=" * 60)
print("SIMPLE GENERIC KINEMATICS")
print("=" * 60)

def mandelstams_from_dual_momenta(n, seed):
    """
    Use dual momenta in high-dimensional space to get generic s_ij.
    In D >> n dimensions, we can place n null momenta generically.
    
    Method: Use D=n-1 massless momenta with all components random.
    """
    set_random_seed(seed)
    
    D = 10  # High enough dimension
    
    # Generate n-1 random momenta in R^D (not null, just random)
    momenta = []
    for i in range(n-1):
        p = vector(QQ, [randint(-10, 10) for _ in range(D)])
        momenta.append(p)
    
    # Last momentum from conservation
    momenta.append(-sum(momenta))
    
    # Compute s_ij = 2 p_i . p_j with Euclidean metric (or Minkowski)
    # For genericity, use Euclidean
    s = {}
    for i in range(n):
        for j in range(i+1, n):
            s[(i+1, j+1)] = 2 * momenta[i].dot_product(momenta[j])
    
    # Check momentum conservation
    for i in range(1, n+1):
        total = sum(s.get((min(i,j), max(i,j)), 0) for j in range(1, n+1) if j != i)
        if total != 0:
            print(f"  WARNING: Particle {i} violation: {total}")
    
    return s


def count_solutions(s, n, F, gauge):
    sigma1, sigma2, sigma3 = gauge
    
    var_names = ['x%d' % i for i in range(4, n+1)]
    R = PolynomialRing(F, var_names)
    xs = R.gens()
    
    sigmas = {1: F(sigma1), 2: F(sigma2), 3: F(sigma3)}
    for i, x in enumerate(xs):
        sigmas[i+4] = x
    
    def get_s(i, j):
        return F(s[(min(i,j), max(i,j))])
    
    polys = []
    for a in range(4, n+1):
        h_cleared = R(0)
        for b in range(1, n+1):
            if b != a:
                term = get_s(a, b)
                for c in range(1, n+1):
                    if c != a and c != b:
                        term *= (sigmas[a] - sigmas[c])
                h_cleared += term
        polys.append(h_cleared)
    
    I = R.ideal(polys)
    
    S = R(1)
    for x in xs:
        S *= x * (x - F(sigma1)) * (x - F(sigma2)) * (x - F(sigma3))
    for i in range(len(xs)):
        for j in range(i+1, len(xs)):
            S *= (xs[i] - xs[j])
    
    I_sat, _ = I.saturation(S)
    sols = I_sat.variety()
    return len(sols)


n = 7
gauge = (0, 1, 2)
seed = int(sys.argv[1]) if len(sys.argv) > 1 else 0
p = int(sys.argv[2]) if len(sys.argv) > 2 else 101

print(f"n = {n}, seed = {seed}, p = {p}")

s = mandelstams_from_dual_momenta(n, seed)
print("Mandelstams generated")

# Check non-degeneracy: no s_ij should be zero
zeros = [(i,j) for (i,j), v in s.items() if v == 0]
if zeros:
    print(f"  Zero mandelstams: {zeros}")
else:
    print("  No zero mandelstams - GOOD")

print(f"\nCounting over F_{p}...")
t0 = time.time()
count = count_solutions(s, n, GF(p), gauge)
t1 = time.time()

print(f"\n*** {count} SOLUTIONS ({t1-t0:.1f}s) ***")

if count == 24:
    print("\n***** 24 SOLUTIONS - EXACTLY (n-3)! *****")
elif count > 20:
    print(f"\n*** CLOSE TO 24! ({24-count} missing) ***")
