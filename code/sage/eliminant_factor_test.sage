#!/usr/bin/env sage
"""
Factor the eliminant using working code from verify_n7_np2_extension.sage
"""

import sys
import time
from datetime import datetime

# =============================================================================
# CLI
# =============================================================================

def _get_arg(flag, default=None):
    if flag in sys.argv:
        i = sys.argv.index(flag)
        if i + 1 < len(sys.argv):
            return sys.argv[i + 1]
    return default

seed = int(_get_arg("--seed", "0"))
p = int(_get_arg("--p", "17"))

print("=" * 60)
print("ELIMINANT FACTORIZATION TEST")
print(f"seed={seed}, p={p}")
print("=" * 60)

# =============================================================================
# KINEMATICS (from verify_n7_np2_extension.sage)
# =============================================================================

def generate_4D_null_momenta(n, seed=42, scale=100):
    set_random_seed(seed)
    
    def random_null_vector(scale):
        while True:
            m = randint(1, scale)
            n_val = randint(1, scale)
            p_var = randint(1, scale)
            q = randint(1, scale)
            
            x = m^2 + n_val^2 - p_var^2 - q^2
            y = 2*(m*q + n_val*p_var)
            z = 2*(n_val*q - m*p_var)
            E = m^2 + n_val^2 + p_var^2 + q^2
            
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
            
            def is_null(pvec):
                return pvec[0]^2 - pvec[1]^2 - pvec[2]^2 - pvec[3]^2 == 0
            
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


def good_reduction(mandelstams, p):
    for key, val in mandelstams.items():
        if hasattr(val, "denominator") and (val.denominator() % p == 0):
            return False
    return True


# =============================================================================
# SOLUTION FINDING WITH ELIMINANT
# =============================================================================

def find_solutions_and_eliminant(mandelstams, gauge, F):
    sigma1, sigma2, sigma3 = gauge
    
    R = PolynomialRing(F, 'x4,x5,x6,x7')
    x4, x5, x6, x7 = R.gens()
    
    sigmas = {1: F(sigma1), 2: F(sigma2), 3: F(sigma3), 4: x4, 5: x5, 6: x6, 7: x7}
    
    def get_s(i, j):
        return F(mandelstams[(min(i,j), max(i,j))])
    
    polys = []
    for a in range(4, 8):
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
    
    # Spurious factors
    S = x4 * (x4 - F(sigma1)) * (x4 - F(sigma2)) * (x4 - F(sigma3))
    S *= x5 * (x5 - F(sigma1)) * (x5 - F(sigma2)) * (x5 - F(sigma3))
    S *= x6 * (x6 - F(sigma1)) * (x6 - F(sigma2)) * (x6 - F(sigma3))
    S *= x7 * (x7 - F(sigma1)) * (x7 - F(sigma2)) * (x7 - F(sigma3))
    S *= (x4 - x5) * (x4 - x6) * (x4 - x7)
    S *= (x5 - x6) * (x5 - x7)
    S *= (x6 - x7)
    
    print("  Saturating ideal...")
    t0 = time.time()
    I_sat, _ = I.saturation(S)
    t1 = time.time()
    print(f"  Saturation done in {t1-t0:.1f}s")
    
    print("  Computing variety...")
    t0 = time.time()
    sols = I_sat.variety()
    t1 = time.time()
    print(f"  Found {len(sols)} solutions in {t1-t0:.1f}s")
    
    if len(sols) == 0:
        return [], None
    
    # Extract x4 values (eliminant in x4)
    x4_vals = [sol[x4] for sol in sols]
    
    # Build eliminant
    Rx = PolynomialRing(F, 'X')
    X = Rx.gen()
    elim = prod(X - v for v in x4_vals)
    
    return sols, elim


# =============================================================================
# MAIN
# =============================================================================

n = 7
print(f"\nGenerating kinematics for n={n}, seed={seed}...")
momenta = generate_4D_null_momenta(n, seed)
mandelstams = momenta_to_mandelstams(momenta)

# Check good reduction
if not good_reduction(mandelstams, p):
    print(f"Bad reduction at p={p}, trying p+4...")
    p = p + 4
    while not good_reduction(mandelstams, p):
        p += 4

print(f"Using prime p={p}")

# Standard gauge
gauge = (0, 1, 2)  # sigma_1=0, sigma_2=1, sigma_3=2

F = GF(p)
print(f"\nFinding solutions over F_{p}...")

sols, elim = find_solutions_and_eliminant(mandelstams, gauge, F)

if elim is not None:
    print(f"\n" + "=" * 60)
    print("ELIMINANT ANALYSIS")
    print("=" * 60)
    print(f"Degree: {elim.degree()}")
    
    print(f"\nFactoring...")
    fact = elim.factor()
    print(f"Factorization: {fact}")
    
    degrees = []
    for f, e in fact:
        degrees.extend([f.degree()] * e)
    degrees.sort(reverse=True)
    
    print(f"\nFactor degrees: {degrees}")
    print(f"Partition: {' + '.join(map(str, degrees))} = {sum(degrees)}")
    
    # Analysis
    print("\n" + "=" * 60)
    print("INTERPRETATION")
    print("=" * 60)
    
    if sum(degrees) == 24:
        if degrees.count(12) == 2 and len(degrees) == 2:
            print("*** 12 + 12 FACTORIZATION! Two degree-12 factors ***")
            print("This supports the hypothesis of two Galois orbits with group 12T48")
        elif max(degrees) == 24 and len(degrees) == 1:
            print("*** IRREDUCIBLE! Single degree-24 factor ***")
            print("The Galois group acts transitively on all 24 solutions")
        else:
            print(f"*** INTERESTING FACTORIZATION ***")
            print(f"Factor structure suggests multiple Galois orbits")
    else:
        print(f"*** UNEXPECTED: Sum of degrees is {sum(degrees)}, not 24 ***")
else:
    print("No solutions found over F_p")
    print("Try a different seed or prime")
