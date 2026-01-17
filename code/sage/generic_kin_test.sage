#!/usr/bin/env sage
"""
Test with GENERIC Mandelstam invariants (not 4D-constrained)
This should give the full 24 solutions.
"""

import sys
import time

print("=" * 60)
print("GENERIC KINEMATICS TEST (not 4D-restricted)")
print("=" * 60)

def generate_generic_mandelstams(n, seed):
    """
    Generate generic Mandelstam invariants satisfying only:
      sum_j s_ij = 0 for each i (momentum conservation)
    
    This is what CHY expects - NO 4D Gram constraints.
    
    Method: Generate (n-1)(n-2)/2 - 1 free parameters, then solve for rest.
    Actually, simpler: generate all s_ij freely, then project to constraint surface.
    
    Even simpler: use parameterization that automatically satisfies constraints.
    For n particles, momentum conservation gives n constraints, but they sum to zero,
    so effectively n-1 independent constraints.
    
    Degrees of freedom: C(n,2) - (n-1) = n(n-1)/2 - (n-1) = (n-1)(n-2)/2
    For n=7: 6*5/2 = 15 free parameters.
    """
    set_random_seed(seed)
    
    # Start with all s_ij = 0
    s = {}
    for i in range(1, n+1):
        for j in range(i+1, n+1):
            s[(i,j)] = QQ(0)
    
    # Add random contributions that preserve momentum conservation
    # Each s_ij enters constraint for particle i and particle j
    # We can freely set s_ij for a spanning tree, then others are determined
    
    # Actually, let's use the explicit parameterization:
    # Generate random s_ij for 1 <= i < j <= n-1, then solve for s_{i,n}
    
    for i in range(1, n):
        for j in range(i+1, n):
            s[(i,j)] = QQ(randint(-50, 50))
            if s[(i,j)] == 0:
                s[(i,j)] = QQ(randint(1, 50))
    
    # Solve for s_{i,n} using sum_j s_ij = 0 for i = 1, ..., n-1
    # For particle i (i < n): s_{i,n} = - sum_{j != i, j != n} s_{ij}
    for i in range(1, n):
        total = QQ(0)
        for j in range(1, n):
            if j != i:
                total += s.get((min(i,j), max(i,j)), QQ(0))
        s[(i,n)] = -total
    
    # Check: momentum conservation for particle n should now be automatic
    # sum_{j < n} s_{j,n} = sum_{j < n} (- sum_{k != j, k != n} s_{jk})
    # This is -(sum over all pairs (j,k) with j,k < n of s_{jk} counted twice) = 0? No...
    # Actually: sum_j s_{nj} = sum_{j<n} s_{j,n} = sum_{j<n} (-sum_{k!=j, k<n} s_{jk})
    # Each s_{jk} appears in both -s_{j,n} and -s_{k,n} terms, so cancels.
    # Wait, that's not right either...
    
    # Let me verify explicitly
    
    return s


def count_chy_solutions(s, n, F, gauge):
    """Count CHY solutions over field F."""
    sigma1, sigma2, sigma3 = gauge
    
    # Variables for sigma_4, ..., sigma_n
    var_names = ['x%d' % i for i in range(4, n+1)]
    R = PolynomialRing(F, var_names)
    xs = R.gens()
    
    sigmas = {1: F(sigma1), 2: F(sigma2), 3: F(sigma3)}
    for i, x in enumerate(xs):
        sigmas[i+4] = x
    
    def get_s(i, j):
        key = (min(i,j), max(i,j))
        val = s.get(key, 0)
        return F(val)
    
    # Build scattering equations (cleared denominators)
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
    
    # Spurious solutions: collisions with gauge points or each other
    S = R(1)
    for x in xs:
        S *= x * (x - F(sigma1)) * (x - F(sigma2)) * (x - F(sigma3))
    for i in range(len(xs)):
        for j in range(i+1, len(xs)):
            S *= (xs[i] - xs[j])
    
    I_sat, _ = I.saturation(S)
    
    try:
        sols = I_sat.variety()
        return len(sols)
    except Exception as e:
        print(f"  Error: {e}")
        return -1


# Parameters
n = 7
gauge = (0, 1, 2)
seed = int(sys.argv[1]) if len(sys.argv) > 1 else 0
p = int(sys.argv[2]) if len(sys.argv) > 2 else 101

print(f"n = {n}, seed = {seed}, p = {p}")
print()

# Generate generic kinematics
print("Generating GENERIC Mandelstams (momentum conservation only)...")
s = generate_generic_mandelstams(n, seed)

# Verify momentum conservation
print("Verifying momentum conservation...")
for i in range(1, n+1):
    total = sum(s.get((min(i,j), max(i,j)), 0) for j in range(1, n+1) if j != i)
    if total != 0:
        print(f"  WARNING: sum_j s_{i}j = {total} != 0")
    else:
        print(f"  Particle {i}: OK")

# Count solutions
print(f"\nCounting solutions over F_{p}...")
t0 = time.time()
F = GF(p)
count = count_chy_solutions(s, n, F, gauge)
t1 = time.time()

print(f"\n*** RESULT: {count} solutions over F_{p} ({t1-t0:.1f}s) ***")

if count == 24:
    print("\n***** EXACTLY 24 SOLUTIONS! *****")
    print("***** GENERIC KINEMATICS WORK! *****")
elif count == (n-3):
    print(f"\n(n-3)! = {factorial(n-3)} expected")
    print(f"Got {count} - this is {count} / 24 = {count/24:.2%}")
else:
    print(f"\nExpected 24, got {count}")
    print("Try different seed or prime")
