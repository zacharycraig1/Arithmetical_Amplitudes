"""
VERIFY LOOP χ₈ CONSTRAINT
=========================

This verifies that the cross-loop product P = A × d_1L × d_2L_s × d_2L_t
satisfies χ_p(P) = χ₈(p) = χ_p(-2).

Uses the explicit BDS formulas from northstar38.
"""

from sage.all import *

def angle_bracket(lam, i, j):
    """Compute angle bracket <ij>"""
    return lam[i][0] * lam[j][1] - lam[i][1] * lam[j][0]

def square_bracket(lam_tilde, i, j):
    """Compute square bracket [ij]"""
    return lam_tilde[i][0] * lam_tilde[j][1] - lam_tilde[i][1] * lam_tilde[j][0]

def mandelstam_s(lam, lam_tilde, i, j):
    """Compute s_ij = <ij>[ji]"""
    return angle_bracket(lam, i, j) * square_bracket(lam_tilde, j, i)

def A4_tree_MHV(lam, F):
    """
    4-point MHV tree amplitude: A = <12>³ / (<23><34><41>)
    """
    ab12 = angle_bracket(lam, 0, 1)
    ab23 = angle_bracket(lam, 1, 2)
    ab34 = angle_bracket(lam, 2, 3)
    ab41 = angle_bracket(lam, 3, 0)
    
    denom = ab23 * ab34 * ab41
    if denom == F(0):
        return None
    
    return ab12**3 / denom

def get_mandelstams(lam, lam_tilde):
    """Get s and t"""
    s = mandelstam_s(lam, lam_tilde, 0, 1)
    t = mandelstam_s(lam, lam_tilde, 1, 2)
    return s, t

def compute_cross_loop_product(lam, lam_tilde, F):
    """
    Compute P = A × d_1L × d_2L_s × d_2L_t
    
    d_1L   = -(1/2) s t A
    d_2L_s = +(1/4) s² t A
    d_2L_t = +(1/4) s t² A
    
    P = A × (-½st·A) × (¼s²t·A) × (¼st²·A)
      = A⁴ × (-½) × (¼)² × s⁴ t⁴
      = A⁴ s⁴ t⁴ × (-1/32)
    """
    A = A4_tree_MHV(lam, F)
    if A is None or A == F(0):
        return None
    
    s, t = get_mandelstams(lam, lam_tilde)
    if s == F(0) or t == F(0):
        return None
    
    # Compute each coefficient
    half = F(2).inverse()
    quarter = F(4).inverse()
    
    d_1L = -half * s * t * A
    d_2L_s = quarter * s**2 * t * A
    d_2L_t = quarter * s * t**2 * A
    
    # The product
    P = A * d_1L * d_2L_s * d_2L_t
    
    return P

def test_chi8_constraint(num_samples=100, seed=42):
    """Test that χ_p(P) = χ₈(p) for various primes and kinematics."""
    
    print("""
================================================================================
             VERIFY LOOP χ₈ CONSTRAINT
================================================================================

Testing: χ_p(P) = χ₈(p) = χ_p(-2)

Where P = A × d_1L × d_2L_s × d_2L_t = (Ast)⁴ / (-32)
""")
    
    primes = [101, 103, 107, 109, 113, 127, 131, 137, 139, 149]
    
    results = []
    
    for p in primes:
        F = GF(p)
        chi8_p = kronecker_symbol(-2, p)
        
        set_random_seed(seed)
        
        valid = 0
        match = 0
        
        for _ in range(num_samples):
            # Random spinors
            lam = [[F.random_element(), F.random_element()] for _ in range(4)]
            lam_tilde = [[F.random_element(), F.random_element()] for _ in range(4)]
            
            P = compute_cross_loop_product(lam, lam_tilde, F)
            
            if P is None or P == F(0):
                continue
            
            valid += 1
            chi_P = kronecker_symbol(int(P), p)
            
            if chi_P == chi8_p:
                match += 1
        
        match_rate = float(match / valid) if valid > 0 else 0.0
        results.append((p, chi8_p, valid, match, match_rate))
    
    print(f"{'Prime':>6} | {'χ₈(p)':>6} | {'Valid':>6} | {'Match':>6} | {'Rate':>8}")
    print("-" * 45)
    
    total_valid = 0
    total_match = 0
    
    for p, chi8_p, valid, match, rate in results:
        print(f"{p:>6} | {chi8_p:>6} | {valid:>6} | {match:>6} | {rate*100:>7.1f}%")
        total_valid += valid
        total_match += match
    
    print("-" * 45)
    overall_rate = float(total_match / total_valid) if total_valid > 0 else 0.0
    print(f"{'TOTAL':>6} | {'':>6} | {total_valid:>6} | {total_match:>6} | {overall_rate*100:>7.1f}%")
    
    print()
    if total_match == total_valid:
        print("*** PERFECT 100% MATCH ***")
        print("*** χ₈ CONSTRAINT VERIFIED ***")
    else:
        violations = total_valid - total_match
        print(f"*** {violations} VIOLATIONS FOUND ***")
    
    return results

def verify_algebraic_identity():
    """Verify the algebraic identity [P] = [-2]."""
    
    print("""
================================================================================
             ALGEBRAIC IDENTITY VERIFICATION
================================================================================

Proving: P = (Ast)⁴ / (-32) ≡ -2 in K*/(K*)²
""")
    
    R = QQ['A', 's', 't']
    A, s, t = R.gens()
    
    # The product
    P_num = A**4 * s**4 * t**4
    P_denom = -32
    
    print(f"P = {P_num} / ({P_denom})")
    print()
    
    # Factor out perfect squares
    # (Ast)⁴ = ((Ast)²)² is a perfect square
    # -32 = -1 × 32 = -1 × 2⁵ = -1 × 2 × 2⁴ = -1 × 2 × (2²)²
    
    print("Factorization:")
    print("  (Ast)⁴ = ((Ast)²)² is a perfect square → [1]")
    print("  -32 = -1 × 2⁵ = -1 × 2 × (2²)²")
    print("      = -1 × 2 × (perfect square)")
    print("      → [-1] × [2] = [-2]")
    print()
    print("Therefore: [P] = [1] × [-2] = [-2]")
    print()
    print("Finite field consequence:")
    print("  χ_p(P) = χ_p(-2) = χ_p(-1) × χ_p(2) = χ₈(p)")

if __name__ == "__main__":
    verify_algebraic_identity()
    print()
    test_chi8_constraint(num_samples=200)
