"""
Kinematics Generator for n-point 4D Null Momenta
=================================================

This module provides deterministic generation of rational 4D null momenta
satisfying momentum conservation for CHY scattering equation verification.

Method: Pythagorean null vector parametrization with constrained final momenta.
"""

def make_null_vector(m, n, p, q, signs=(1, 1, 1)):
    """
    Generate a null 4-vector using the Pythagorean parametrization.
    
    Returns v = (E, x, y, z) satisfying E^2 - x^2 - y^2 - z^2 = 0.
    
    Formula:
        E = m^2 + n^2 + p^2 + q^2
        x = (m^2 + n^2 - p^2 - q^2) * signs[0]
        y = 2*(m*q + n*p) * signs[1]
        z = 2*(n*q - m*p) * signs[2]
    """
    E = m^2 + n^2 + p^2 + q^2
    x = (m^2 + n^2 - p^2 - q^2) * signs[0]
    y = 2*(m*q + n*p) * signs[1]
    z = 2*(n*q - m*p) * signs[2]
    return (QQ(E), QQ(x), QQ(y), QQ(z))


def check_null(p):
    """Verify p^2 = 0 in Minkowski signature (+,-,-,-)."""
    return p[0]^2 - p[1]^2 - p[2]^2 - p[3]^2 == 0


def check_momentum_conservation(momenta):
    """Verify sum of all momenta is zero."""
    n = len(momenta)
    for mu in range(4):
        total = sum(momenta[a][mu] for a in range(n))
        if total != 0:
            return False
    return True


def levi_civita(p1, p2, p3, p4):
    """Compute epsilon(p1, p2, p3, p4) = det([p1; p2; p3; p4])."""
    M = matrix(QQ, [list(p1), list(p2), list(p3), list(p4)])
    return M.det()


def make_4d_massless_point(n, seed, scale=100):
    """
    Generate n null 4D momenta with momentum conservation.
    
    Algorithm:
    1. Generate n-2 random null vectors using Pythagorean parametrization
    2. Solve for last 2 momenta to enforce sum = 0
    3. Verify all constraints
    
    Returns: list of n tuples (E, x, y, z)
    """
    set_random_seed(seed)
    
    # Generate n-2 random null momenta
    momenta = []
    for _ in range(n - 2):
        # Random integers for Pythagorean parametrization
        m = randint(1, scale)
        n_val = randint(1, scale)
        p = randint(1, scale)
        q = randint(1, scale)
        signs = (choice([-1, 1]), choice([-1, 1]), choice([-1, 1]))
        v = make_null_vector(m, n_val, p, q, signs)
        momenta.append(v)
    
    # Sum of first n-2 momenta
    K = [sum(momenta[a][mu] for a in range(n-2)) for mu in range(4)]
    K_sq = K[0]^2 - K[1]^2 - K[2]^2 - K[3]^2
    
    # If K is null, we need to regenerate
    if K_sq == 0:
        return make_4d_massless_point(n, seed + 1000, scale)
    
    # Find p_{n-1} and p_n such that:
    # - p_{n-1}^2 = 0
    # - p_n^2 = 0
    # - p_{n-1} + p_n = -K
    
    # We parameterize: p_{n-1} has components (E, x, y, z) with E^2 = x^2 + y^2 + z^2
    # and p_n = -K - p_{n-1} must also be null
    
    for attempt in range(1000):
        # Random y, z for p_{n-1}
        y = QQ(randint(-scale, scale))
        z = QQ(randint(-scale, scale))
        
        if K[0] == 0:
            continue
        
        # Solve for x from the constraint that p_n is null
        a = K[1]^2 - K[0]^2
        if a == 0:
            continue
        
        C = y*K[2] + z*K[3] - K_sq/2
        b = 2*K[1]*C
        c_val = C^2 - K[0]^2*(y^2 + z^2)
        
        disc = b^2 - 4*a*c_val
        if disc >= 0 and disc.is_square():
            x = (-b + disc.sqrt()) / (2*a)
            E = (x*K[1] + y*K[2] + z*K[3] - K_sq/2) / K[0]
            
            p_n_minus_1 = (E, x, y, z)
            p_n = tuple(-K[mu] - p_n_minus_1[mu] for mu in range(4))
            
            if check_null(p_n_minus_1) and check_null(p_n):
                full_momenta = momenta + [p_n_minus_1, p_n]
                
                # Verify all constraints
                assert all(check_null(p) for p in full_momenta), "Masslessness failed"
                assert check_momentum_conservation(full_momenta), "Momentum conservation failed"
                
                # Check genericity
                eps = levi_civita(full_momenta[0], full_momenta[1], 
                                  full_momenta[2], full_momenta[3])
                if eps == 0:
                    continue  # Degenerate, try again
                
                return full_momenta
    
    # Failed to find, try different seed
    return make_4d_massless_point(n, seed + 1000, scale)


def momenta_to_mandelstams(momenta):
    """Compute Mandelstam invariants s_{ij} = 2 p_i · p_j."""
    n = len(momenta)
    s = {}
    for i in range(n):
        for j in range(i+1, n):
            pi, pj = momenta[i], momenta[j]
            dot = pi[0]*pj[0] - pi[1]*pj[1] - pi[2]*pj[2] - pi[3]*pj[3]
            s[(i+1, j+1)] = 2 * dot
    return s


def check_no_decoupling(mandelstams, n, p):
    """
    Check that no particle is decoupled mod p.
    
    A particle a is decoupled if s_{ab} ≡ 0 (mod p) for all b ≠ a.
    """
    F = GF(p)
    for a in range(1, n+1):
        all_zero = True
        for b in range(1, n+1):
            if a != b:
                key = (min(a, b), max(a, b))
                if F(mandelstams[key]) != 0:
                    all_zero = False
                    break
        if all_zero:
            return False, a
    return True, None


def good_reduction(mandelstams, n, p):
    """
    Check if prime p has good reduction for the kinematic point.
    
    Conditions:
    1. No Mandelstam denominator vanishes mod p
    2. No particle is decoupled mod p
    """
    # Check denominators
    for key, val in mandelstams.items():
        if hasattr(val, 'denominator') and val.denominator() % p == 0:
            return False, "denominator"
    
    # Check decoupling
    ok, particle = check_no_decoupling(mandelstams, n, p)
    if not ok:
        return False, f"particle {particle} decoupled"
    
    return True, "good"


# =============================================================================
# SANITY CHECK (only runs when script is executed directly, not load()'ed)
# =============================================================================

def _run_tests():
    print("Testing kinematics generator...")
    
    for seed in [42, 100, 200]:
        momenta = make_4d_massless_point(7, seed)
        s = momenta_to_mandelstams(momenta)
        
        print(f"\nSeed {seed}:")
        print(f"  Momenta generated: {len(momenta)}")
        print(f"  All null: {all(check_null(p) for p in momenta)}")
        print(f"  Conservation: {check_momentum_conservation(momenta)}")
        
        eps = levi_civita(momenta[0], momenta[1], momenta[2], momenta[3])
        print(f"  ε(p1,p2,p3,p4) = {eps} (nonzero: {eps != 0})")
        
        for p in [7, 11, 19, 23]:
            ok, reason = good_reduction(s, 7, p)
            print(f"  p={p}: {reason}")
    
    print("\n✓ Kinematics generator working correctly")


# Only run tests when executed directly, not when load()'ed
import sys
if __name__ == "__main__" and len(sys.argv) > 0 and "kinematics" in sys.argv[0]:
    _run_tests()
