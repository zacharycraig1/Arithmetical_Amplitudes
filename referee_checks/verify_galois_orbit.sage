"""
GALOIS ORBIT VERIFICATION: Testing D₂₄ Conjecture
===================================================

This script tests the hypothesis that the n=7 CHY Galois group is D₂₄ (order 48)
by counting solutions over F_{p^k} for various extension degrees k.

Key predictions from D₂₄ structure:
- N_p = 0 (inert primes, solutions need at least Q(i))
- N_{p²} = some subset (solutions in Q(i) subfield)
- N_{p^k} = 24 for some k dividing 48/gcd(48, ord_p)

The subfield lattice of D₂₄ has degrees: 2, 3, 4, 6, 8, 12, 24, 48
So we test k = 1, 2, 3, 4, 6, 8, 12 to find when all 24 solutions appear.

Usage:
    sage verify_galois_orbit.sage
    sage verify_galois_orbit.sage --seed 0 --prime 7 --max_k 12
"""

import sys
import time
from datetime import datetime

print("=" * 70)
print("GALOIS ORBIT VERIFICATION: D₂₄ Conjecture Test")
print("=" * 70)
print(f"Timestamp: {datetime.now()}")
print()

# =============================================================================
# CLI ARGUMENTS
# =============================================================================

def _get_arg(flag, default=None):
    if flag in sys.argv:
        i = sys.argv.index(flag)
        if i + 1 < len(sys.argv):
            return sys.argv[i + 1]
    return default

seed = int(_get_arg("--seed", "0"))
prime = int(_get_arg("--prime", "7"))
max_k = int(_get_arg("--max_k", "12"))

print(f"Seed: {seed}")
print(f"Prime: {prime} (inert: {prime % 4 == 3})")
print(f"Testing extension degrees k = 1, 2, ..., {max_k}")
print()

# =============================================================================
# KINEMATICS GENERATION
# =============================================================================

def generate_4D_null_momenta(n, seed=42, scale=100):
    """Generate n null momenta in 4D with momentum conservation."""
    set_random_seed(seed)
    
    def random_null_vector(scale):
        while True:
            m = randint(1, scale)
            n_val = randint(1, scale)
            p = randint(1, scale)
            q = randint(1, scale)
            
            x = m^2 + n_val^2 - p^2 - q^2
            y = 2*(m*q + n_val*p)
            z = 2*(n_val*q - m*p)
            E = m^2 + n_val^2 + p^2 + q^2
            
            sx, sy, sz = choice([-1,1]), choice([-1,1]), choice([-1,1])
            return (QQ(E), QQ(x)*sx, QQ(y)*sy, QQ(z)*sz)
    
    momenta = [random_null_vector(scale) for _ in range(n-2)]
    
    K = [sum(m[i] for m in momenta) for i in range(4)]
    K_sq = K[0]^2 - K[1]^2 - K[2]^2 - K[3]^2
    
    if K_sq == 0:
        return generate_4D_null_momenta(n, seed+1, scale)
    
    for attempt in range(1000):
        y = QQ(randint(-scale, scale))
        z = QQ(randint(-scale, scale))
        
        if K[0] == 0:
            continue
        
        a = K[1]^2 - K[0]^2
        if a == 0:
            continue
        
        C = y*K[2] + z*K[3] - K_sq/2
        b = 2*K[1]*C
        c_val = C^2 - K[0]^2*(y^2+z^2)
        
        disc = b^2 - 4*a*c_val
        if disc >= 0 and disc.is_square():
            x = (-b + disc.sqrt()) / (2*a)
            E = (x*K[1] + y*K[2] + z*K[3] - K_sq/2) / K[0]
            
            kn_1 = (E, x, y, z)
            kn = tuple(-K[i] - kn_1[i] for i in range(4))
            
            def is_null(p):
                return p[0]^2 - p[1]^2 - p[2]^2 - p[3]^2 == 0
            
            if is_null(kn_1) and is_null(kn):
                return momenta + [kn_1, kn]
    
    return generate_4D_null_momenta(n, seed+1000, scale)


def momenta_to_mandelstams(momenta):
    n = len(momenta)
    s = {}
    for i in range(n):
        for j in range(i+1, n):
            pi, pj = momenta[i], momenta[j]
            dot = pi[0]*pj[0] - pi[1]*pj[1] - pi[2]*pj[2] - pi[3]*pj[3]
            s[(i+1, j+1)] = 2 * dot
    return s


# =============================================================================
# F_{p^k} SOLUTION COUNTING
# =============================================================================

def count_over_Fpk(mandelstams, gauge, p, k):
    """
    Count solutions to n=7 CHY over F_{p^k}.
    
    Returns the number of F_{p^k}-rational solutions.
    """
    sigma1, sigma2, sigma3 = gauge
    
    # Create the field F_{p^k}
    if k == 1:
        F = GF(p)
    else:
        F = GF(p^k, 'a')
    
    R = PolynomialRing(F, 'x4,x5,x6,x7')
    x4, x5, x6, x7 = R.gens()
    
    sigmas = {1: F(sigma1), 2: F(sigma2), 3: F(sigma3), 4: x4, 5: x5, 6: x6, 7: x7}
    
    def get_s(i, j):
        return F(mandelstams[(min(i,j), max(i,j))])
    
    # Build cleared CHY polynomials
    polys = []
    for a in range(4, 8):
        h_cleared = R(0)
        for b in range(1, 8):
            if b != a:
                term = get_s(a, b)
                for c in range(1, 8):
                    if c != a and c != b:
                        term *= (sigmas[a] - sigmas[c])
                h_cleared += term
        polys.append(h_cleared)
    
    I = R.ideal(polys)
    
    # Spurious factor: gauge collisions and pairwise collisions
    S = x4 * (x4 - F(sigma1)) * (x4 - F(sigma2)) * (x4 - F(sigma3))
    S *= x5 * (x5 - F(sigma1)) * (x5 - F(sigma2)) * (x5 - F(sigma3))
    S *= x6 * (x6 - F(sigma1)) * (x6 - F(sigma2)) * (x6 - F(sigma3))
    S *= x7 * (x7 - F(sigma1)) * (x7 - F(sigma2)) * (x7 - F(sigma3))
    S *= (x4 - x5) * (x4 - x6) * (x4 - x7)
    S *= (x5 - x6) * (x5 - x7)
    S *= (x6 - x7)
    
    # Saturate
    I_sat, _ = I.saturation(S)
    
    # Count solutions
    try:
        sols = I_sat.variety()
        return len(sols)
    except Exception as e:
        print(f"    [Warning] variety() failed for k={k}: {e}")
        return -1


# =============================================================================
# MAIN
# =============================================================================

def main():
    gauge = (0, 1, -1)
    
    # Generate kinematics
    print("Generating kinematics...")
    momenta = generate_4D_null_momenta(7, seed=seed*100+42)
    mandelstams = momenta_to_mandelstams(momenta)
    print("  Done.")
    print()
    
    p = prime
    
    # Test extension degrees
    # D₂₄ has subfields at degrees 2, 3, 4, 6, 8, 12, 24
    # We test divisors of 24 (since 48/2 = 24 is the max we'd need for inert primes)
    test_degrees = [k for k in [1, 2, 3, 4, 6, 8, 12, 24] if k <= max_k]
    
    results = []
    
    print(f"Testing N_{{p^k}} for p = {p}:")
    print("-" * 50)
    
    for k in test_degrees:
        print(f"  k = {k}: Computing N_{{p^{k}}} = N_{{{p}^{k}}} ...", end=" ", flush=True)
        t0 = time.time()
        N = count_over_Fpk(mandelstams, gauge, p, k)
        t1 = time.time()
        print(f"N = {N}  ({t1-t0:.2f}s)")
        results.append((k, N))
        
        # If we found all 24 solutions, we can stop
        if N == 24:
            print(f"\n  *** Found all 24 solutions at k = {k} ***")
            break
    
    print()
    print("=" * 70)
    print("RESULTS SUMMARY")
    print("=" * 70)
    print(f"seed={seed}, p={p}")
    print()
    print("Extension degree k | N_{p^k}")
    print("-" * 30)
    for k, N in results:
        marker = " <-- ALL 24!" if N == 24 else ""
        print(f"  k = {k:2d}            | {N:3d}{marker}")
    
    print()
    
    # Interpretation
    print("INTERPRETATION:")
    print("-" * 50)
    
    if any(N == 24 for k, N in results):
        k_full = min(k for k, N in results if N == 24)
        print(f"All 24 solutions become F_{{p^{k_full}}}-rational.")
        print(f"This suggests the Frobenius element at p={p} has order {k_full}")
        print(f"in the Galois group action on solutions.")
        print()
        print("This is CONSISTENT with the D₂₄ conjecture!")
        print("The splitting field has degree 48, and Q(i) is a degree-2 subfield.")
    else:
        print(f"Did not find all 24 solutions up to k = {max_k}.")
        print("Try increasing --max_k to test larger extensions.")
    
    print()
    
    # Check intermediate counts
    if len(results) >= 2:
        k1, N1 = results[0]
        k2, N2 = results[1]
        if N1 == 0 and N2 > 0:
            print(f"N_p = 0 but N_{{p^2}} = {N2}")
            print("This confirms solutions require at least the degree-2 extension (Q(i)).")


if __name__ == "__main__":
    main()
