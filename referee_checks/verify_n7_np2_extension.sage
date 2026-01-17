"""
REFEREE VERIFICATION: n=7 Quadratic Extension Check
=============================================================

This script verifies the evidence for a Q(i) obstruction:
    - At inert prime p: N_p = 0 (solutions vanish over F_p)
    - Over extension: N_{p^2} > 0 (solutions reappear over F_{p^2})

Interpretation:
- If N_p = 0 but N_{p^2} > 0, this is consistent with a quadratic obstruction
  (Q(i) as a subfield mechanism).
- If additionally N_{p^2} = 24 for many good primes, that would support the stronger
  (and typically false here) statement that the splitting field is exactly Q(i).

Usage:
    sage verify_n7_np2_extension.sage
    sage verify_n7_np2_extension.sage --seed 0 --prime 7
    sage verify_n7_np2_extension.sage --seed_value 42 --prime 7

Expected output:
    N_p = 0, N_{p^2} > 0

WARNING: This is computationally expensive. Use small primes (7, 11) only.
"""

import sys
import time
import json
from datetime import datetime

load("code/sage/kinematics.sage")

print("=" * 70)
print("REFEREE VERIFICATION: n=7 Quadratic Extension Check")
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

seed_id = int(_get_arg("--seed", "0"))
seed_value_override = _get_arg("--seed_value", None)
prime = int(_get_arg("--prime", "7"))

def load_d1_seed(seed_id, seed_value_override=None):
    if seed_value_override is not None:
        return int(seed_value_override)
    with open("data/D1_seeds.json", "r") as handle:
        data = json.load(handle)
    for entry in data.get("seeds", []):
        if int(entry.get("id", -1)) == seed_id:
            return int(entry["seed_value"])
    raise ValueError(f"Seed id {seed_id} not found in data/D1_seeds.json")

seed_value = load_d1_seed(seed_id, seed_value_override)

print(f"Seed id: {seed_id}")
print(f"Seed value (D1): {seed_value}")
print(f"Prime: {prime} (inert: {prime % 4 == 3})")
print()

def good_reduction_or_die(mandelstams, p):
    ok, reason = good_reduction(mandelstams, 7, p)
    if not ok:
        raise ValueError(f"Bad reduction: {reason} at p={p}")


# =============================================================================
# F_q SOLUTION COUNTING (Groebner method)
# =============================================================================

def count_over_Fq(mandelstams, gauge, F, do_variety=True):
    """
    Count solutions to n=7 CHY over field F using Groebner basis.
    
    F can be GF(p) or GF(p^2).
    """
    sigma1, sigma2, sigma3 = gauge
    
    R = PolynomialRing(F, 'x4,x5,x6,x7')
    x4, x5, x6, x7 = R.gens()
    
    sigmas = {1: F(sigma1), 2: F(sigma2), 3: F(sigma3), 4: x4, 5: x5, 6: x6, 7: x7}
    
    def get_s(i, j):
        return F(mandelstams[(min(i,j), max(i,j))])
    
    # Build cleared CHY polynomials
    polys = []
    for a in range(4, 8):
        # h_a = sum_{b != a} s_{ab} / (sigma_a - sigma_b) = 0
        # Cleared: prod_{b != a} (sigma_a - sigma_b) * h_a = 0
        
        denom_product = R(1)
        for b in range(1, 8):
            if b != a:
                denom_product *= (sigmas[a] - sigmas[b])
        
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

    # Optional diagnostic: quotient dimension (may not be available in all Sage builds/types)
    try:
        Q = R.quotient(I_sat)
        dimA = int(Q.vector_space_dimension())
        print(f"  quotient_dim = {dimA}")
    except Exception as e:
        print(f"  Quotient dimension unavailable: {e}")

    # Enumerate solutions (correct point count)
    try:
        sols = I_sat.variety()
        return len(sols)
    except Exception as e:
        print(f"  Variety computation failed: {e}")
        return -1


# =============================================================================
# MAIN
# =============================================================================

def main():
    gauge = (0, 1, -1)
    
    # Generate kinematics (D1 seed)
    print("Generating kinematics (D1)...")
    momenta = make_4d_massless_point(7, seed_value)
    mandelstams = momenta_to_mandelstams(momenta)
    good_reduction_or_die(mandelstams, prime)
    print("  Done.")
    print()
    
    p = prime
    
    # Count over F_p
    print(f"Counting over F_{p}...")
    Fp = GF(p)
    t0 = time.time()
    Np = count_over_Fq(mandelstams, gauge, Fp)
    t1 = time.time()
    print(f"  N_p = {Np}  ({t1-t0:.2f}s)")
    print()
    
    # Count over F_{p^2}
    print(f"Counting over F_{p^2} (this may take a while)...")
    Fp2 = GF(p^2, 'i')
    t2 = time.time()
    Np2 = count_over_Fq(mandelstams, gauge, Fp2, do_variety=True)
    t3 = time.time()
    print(f"  N_{{p^2}} = {Np2}  ({t3-t2:.2f}s)")
    print()
    
    # Summary
    print("=" * 70)
    print("RESULT")
    print("=" * 70)
    print(f"seed_id={seed_id}, seed_value={seed_value}, p={p}:")
    print(f"  N_p      = {Np}   (expected: 0)")
    print(f"  N_{{p^2}} = {Np2}  (expected: >0)")
    print()
    
    if Np == 0 and Np2 > 0:
        print("[PASS] QUADRATIC-EXTENSION WITNESS")
        print("  Solutions vanish over F_p but reappear over F_{p^2}.")
        print("  This is consistent with a Q(i) subfield obstruction mechanism.")
        if Np2 != 24:
            print("  Note: N_{p^2} != 24 indicates the full splitting field is larger than Q(i).")
    else:
        print("[FAIL] Unexpected result")
        sys.exit(1)


if __name__ == "__main__":
    main()
