#!/usr/bin/env sage
"""
Factor eliminant over field extensions to see factorization pattern
"""

import sys
import time

def _get_arg(flag, default=None):
    if flag in sys.argv:
        i = sys.argv.index(flag)
        if i + 1 < len(sys.argv):
            return sys.argv[i + 1]
    return default

seed = int(_get_arg("--seed", "0"))
p = int(_get_arg("--p", "17"))
k = int(_get_arg("--k", "4"))  # Field extension degree

print("=" * 60)
print(f"ELIMINANT OVER F_{{p^{k}}} = F_{{{p}^{k}}}")
print(f"seed={seed}")
print("=" * 60)

# Kinematics
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


def find_solutions(mandelstams, gauge, F):
    sigma1, sigma2, sigma3 = gauge
    
    R = PolynomialRing(F, 'x4,x5,x6,x7')
    x4, x5, x6, x7 = R.gens()
    
    sigmas = {1: F(sigma1), 2: F(sigma2), 3: F(sigma3), 4: x4, 5: x5, 6: x6, 7: x7}
    
    def get_s(i, j):
        return F(mandelstams[(min(i,j), max(i,j))])
    
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
    
    S = x4 * (x4 - F(sigma1)) * (x4 - F(sigma2)) * (x4 - F(sigma3))
    S *= x5 * (x5 - F(sigma1)) * (x5 - F(sigma2)) * (x5 - F(sigma3))
    S *= x6 * (x6 - F(sigma1)) * (x6 - F(sigma2)) * (x6 - F(sigma3))
    S *= x7 * (x7 - F(sigma1)) * (x7 - F(sigma2)) * (x7 - F(sigma3))
    S *= (x4 - x5) * (x4 - x6) * (x4 - x7)
    S *= (x5 - x6) * (x5 - x7)
    S *= (x6 - x7)
    
    I_sat, _ = I.saturation(S)
    sols = I_sat.variety()
    
    return sols, x4


# Main
n = 7
print(f"Generating kinematics...")
momenta = generate_4D_null_momenta(n, seed)
mandelstams = momenta_to_mandelstams(momenta)

gauge = (0, 1, 2)

# Try extension field
print(f"\nComputing over F_{p^k}...")
t0 = time.time()
F = GF(p^k, 'a')
sols, x4_var = find_solutions(mandelstams, gauge, F)
t1 = time.time()

print(f"Found {len(sols)} solutions in {t1-t0:.1f}s")

if len(sols) > 0:
    x4_vals = [sol[x4_var] for sol in sols]
    
    # Minimal polynomial of each x4 value over F_p
    print(f"\nAnalyzing solution structure...")
    
    Fp = GF(p)
    Rx = PolynomialRing(Fp, 'X')
    X = Rx.gen()
    
    # Group by conjugacy classes
    processed = set()
    conjugacy_classes = []
    
    for v in x4_vals:
        if v in processed:
            continue
        
        # Find all F_p conjugates (Frobenius orbit)
        orbit = [v]
        current = v^p
        while current != v:
            orbit.append(current)
            current = current^p
        
        for o in orbit:
            processed.add(o)
        
        conjugacy_classes.append(len(orbit))
    
    conjugacy_classes.sort(reverse=True)
    
    print(f"\nFrobenius orbit sizes (over F_p): {conjugacy_classes}")
    print(f"Total: {sum(conjugacy_classes)}")
    
    # This tells us the degrees of the irreducible factors of the eliminant over F_p
    print(f"\n*** Factor degrees of eliminant over F_p: {conjugacy_classes} ***")
    print(f"*** Partition: {' + '.join(map(str, conjugacy_classes))} = {sum(conjugacy_classes)} ***")
    
    # Check for 12+12 pattern
    if conjugacy_classes == [12, 12]:
        print("\n***** 12 + 12 CONFIRMED! *****")
        print("The eliminant factors into two degree-12 polynomials over Q!")
    elif sum(conjugacy_classes) == 24 and max(conjugacy_classes) == 24:
        print("\n***** IRREDUCIBLE degree-24 *****")
    else:
        print(f"\n***** INTERESTING PATTERN: {conjugacy_classes} *****")
