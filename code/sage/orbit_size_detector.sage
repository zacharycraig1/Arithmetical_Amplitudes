#!/usr/bin/env sage
"""
orbit_size_detector.sage

Compute N_{p^k} for k = 1, 2, 3, 4, 6 at a single inert prime.
This reveals the Frobenius orbit structure directly.

For D24: reflections have sigma^2 = id, so N_{p^2} = 24.
Our data shows N_{p^2} = 4, ruling out D24.

This script answers: What are the orbit sizes?

Usage:
    sage orbit_size_detector.sage --p 7 --seed 0 --max_k 4
"""
import sys
import time
import json
from datetime import datetime

print("=" * 70)
print("ORBIT SIZE DETECTOR: Frobenius orbit analysis over GF(p^k)")
print("=" * 70)

# =============================================================================
# KINEMATICS (adapted from verify_n7_np2_extension.sage)
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
            
            x = m**2 + n_val**2 - p**2 - q**2
            y = 2*(m*q + n_val*p)
            z = 2*(n_val*q - m*p)
            E = m**2 + n_val**2 + p**2 + q**2
            
            sx, sy, sz = choice([-1,1]), choice([-1,1]), choice([-1,1])
            return (QQ(E), QQ(x)*sx, QQ(y)*sy, QQ(z)*sz)
    
    momenta = [random_null_vector(scale) for _ in range(n-2)]
    
    K = [sum(m[i] for m in momenta) for i in range(4)]
    K_sq = K[0]**2 - K[1]**2 - K[2]**2 - K[3]**2
    
    if K_sq == 0:
        return generate_4D_null_momenta(n, seed+1, scale)
    
    for attempt in range(1000):
        y = QQ(randint(-scale, scale))
        z = QQ(randint(-scale, scale))
        
        if K[0] == 0:
            continue
        
        a = K[1]**2 - K[0]**2
        if a == 0:
            continue
        
        C = y*K[2] + z*K[3] - K_sq/2
        b = 2*K[1]*C
        c_val = C**2 - K[0]**2*(y**2+z**2)
        
        disc = b**2 - 4*a*c_val
        if disc >= 0 and disc.is_square():
            x = (-b + disc.sqrt()) / (2*a)
            E = (x*K[1] + y*K[2] + z*K[3] - K_sq/2) / K[0]
            
            kn_1 = (E, x, y, z)
            kn = tuple(-K[i] - kn_1[i] for i in range(4))
            
            def is_null(p):
                return p[0]**2 - p[1]**2 - p[2]**2 - p[3]**2 == 0
            
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
# SOLUTION COUNTING
# =============================================================================

def count_over_Fq(mandelstams, gauge, F, verbose=True):
    """
    Count solutions to n=7 CHY over field F.
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
        denom_product = R(1)
        for b_idx in range(1, 8):
            if b_idx != a:
                denom_product *= (sigmas[a] - sigmas[b_idx])
        
        h_cleared = R(0)
        for b_idx in range(1, 8):
            if b_idx != a:
                term = get_s(a, b_idx)
                for c in range(1, 8):
                    if c != a and c != b_idx:
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
    
    # Enumerate solutions
    try:
        sols = I_sat.variety()
        return len(sols)
    except Exception as e:
        if verbose:
            print(f"    Variety computation failed: {e}")
        return -1


# =============================================================================
# MAIN
# =============================================================================

def _get_arg(flag, default=None):
    if flag in sys.argv:
        i = sys.argv.index(flag)
        if i + 1 < len(sys.argv):
            return sys.argv[i + 1]
    return default


def main():
    seed = int(_get_arg("--seed", "0"))
    p = int(_get_arg("--p", "7"))
    max_k = int(_get_arg("--max_k", "4"))
    
    print(f"\nParameters:")
    print(f"  Prime p = {p}")
    print(f"  Seed = {seed}")
    print(f"  Max extension degree k = {max_k}")
    
    # Check if p is inert
    if p % 4 == 3:
        print(f"  p = 3 mod 4: INERT in Q(i)")
    else:
        print(f"  p = 1 mod 4: SPLIT in Q(i)")
    
    print(f"\nTimestamp: {datetime.now()}")
    print()
    
    # Generate kinematics
    gauge = (0, 1, -1)
    print("Generating kinematics...")
    momenta = generate_4D_null_momenta(7, seed=seed*100+42)
    mandelstams = momenta_to_mandelstams(momenta)
    print("  Done.")
    print()
    
    # Compute N_{p^k} for each k
    results = {}
    
    for k in range(1, max_k + 1):
        print(f"k = {k}: Counting over GF({p}^{k}) = GF({p**k})...")
        F = GF(p**k, 'a')
        t0 = time.time()
        count = count_over_Fq(mandelstams, gauge, F)
        t1 = time.time()
        results[k] = count
        print(f"  N_{{p^{k}}} = {count}  ({t1-t0:.2f}s)")
        print()
    
    # Summary
    print("=" * 70)
    print("SUMMARY: N_{p^k} for k = 1, ..., %d" % max_k)
    print("=" * 70)
    print(f"{'k':>3} | {'p^k':>10} | N_{{p^k}}")
    print("-" * 30)
    for k in range(1, max_k + 1):
        print(f"{k:>3} | {p**k:>10} | {results[k]}")
    
    # Interpretation
    print()
    print("=" * 70)
    print("INTERPRETATION")
    print("=" * 70)
    
    # Find first k where N_{p^k} = 24
    reached_24 = None
    for k in range(1, max_k + 1):
        if results[k] == 24:
            reached_24 = k
            break
    
    if reached_24:
        print(f"All 24 solutions visible over GF(p^{reached_24})")
        print(f"=> Frobenius orbits have sizes dividing {reached_24}")
        
        # Deduce orbit structure from the N_{p^k} sequence
        print()
        print("Orbit structure deduction:")
        for k in range(1, reached_24):
            delta = results.get(k+1, 0) - results.get(k, 0)
            if delta > 0:
                print(f"  {delta} new solutions appear at k={k+1}")
    else:
        print(f"Did not reach 24 solutions by k = {max_k}")
        print("=> Some orbits have size > %d" % max_k)
        print("=> Try larger max_k")
    
    # D24 check
    print()
    print("D24 CHECK:")
    if 2 in results:
        if results[2] == 24:
            print("  N_{p^2} = 24: CONSISTENT with D24 reflections")
        else:
            print(f"  N_{{p^2}} = {results[2]} != 24: RULES OUT D24")
            print("  (D24 reflections have sigma^2 = id, would give N_{p^2} = 24)")
    
    # Save results
    out = {
        "p": int(p),
        "seed": int(seed),
        "max_k": int(max_k),
        "timestamp": str(datetime.now()),
        "results": {int(k): int(v) if isinstance(v, (int, Integer)) else str(v) 
                    for k, v in results.items()}
    }
    outfile = f"orbit_sizes_p{p}_seed{seed}.json"
    with open(outfile, 'w') as f:
        json.dump(out, f, indent=2)
    print(f"\nSaved to {outfile}")


if __name__ == "__main__":
    main()
