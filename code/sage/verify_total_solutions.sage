#!/usr/bin/env sage
"""
Verify the total number of CHY solutions is 24 by computing
the quotient dimension (degree of the 0-dimensional ideal).
"""
import sys
import time

print("=" * 60)
print("CHY n=7: Total Solution Count Verification")
print("=" * 60)

def _get_arg(flag, default=None):
    if flag in sys.argv:
        i = sys.argv.index(flag)
        if i + 1 < len(sys.argv):
            return sys.argv[i + 1]
    return default

seed = int(_get_arg("--seed", "0"))
p = int(_get_arg("--p", "101"))  # Large prime for good reduction

print(f"Prime: {p}")
print(f"Seed: {seed}")
print()

# Generate kinematics
def generate_4D_null_momenta(n, seed=42, scale=100):
    set_random_seed(seed)
    
    def random_null_vector(scale):
        while True:
            m = randint(1, scale)
            n_val = randint(1, scale)
            pv = randint(1, scale)
            q = randint(1, scale)
            
            x = m**2 + n_val**2 - pv**2 - q**2
            y = 2*(m*q + n_val*pv)
            z = 2*(n_val*q - m*pv)
            E = m**2 + n_val**2 + pv**2 + q**2
            
            sx, sy, sz = choice([-1,1]), choice([-1,1]), choice([-1,1])
            return (QQ(E), QQ(x)*sx, QQ(y)*sy, QQ(z)*sz)
    
    momenta = [random_null_vector(scale) for _ in range(n-2)]
    
    K = [sum(mv[i] for mv in momenta) for i in range(4)]
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
            xv = (-b + disc.sqrt()) / (2*a)
            E = (xv*K[1] + y*K[2] + z*K[3] - K_sq/2) / K[0]
            
            kn_1 = (E, xv, y, z)
            kn = tuple(-K[i] - kn_1[i] for i in range(4))
            
            def is_null(pv):
                return pv[0]**2 - pv[1]**2 - pv[2]**2 - pv[3]**2 == 0
            
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


print("Generating kinematics...")
momenta = generate_4D_null_momenta(7, seed=seed*100+42)
mandelstams = momenta_to_mandelstams(momenta)
print("  Done.")

# Build ideal over F_p
F = GF(p)
gauge = (0, 1, -1)
sigma1, sigma2, sigma3 = gauge

R = PolynomialRing(F, 'x4,x5,x6,x7')
x4, x5, x6, x7 = R.gens()

sigmas = {1: F(sigma1), 2: F(sigma2), 3: F(sigma3), 4: x4, 5: x5, 6: x6, 7: x7}

def get_s(i, j):
    return F(mandelstams[(min(i,j), max(i,j))])

# Build CHY polynomials
print("Building CHY ideal...")
polys = []
for a in range(4, 8):
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
print("  Done.")

# Spurious factor
print("Building spurious factor...")
S = x4 * (x4 - F(sigma1)) * (x4 - F(sigma2)) * (x4 - F(sigma3))
S *= x5 * (x5 - F(sigma1)) * (x5 - F(sigma2)) * (x5 - F(sigma3))
S *= x6 * (x6 - F(sigma1)) * (x6 - F(sigma2)) * (x6 - F(sigma3))
S *= x7 * (x7 - F(sigma1)) * (x7 - F(sigma2)) * (x7 - F(sigma3))
S *= (x4 - x5) * (x4 - x6) * (x4 - x7)
S *= (x5 - x6) * (x5 - x7)
S *= (x6 - x7)
print("  Done.")

# Saturate
print("Saturating...")
t0 = time.time()
I_sat, _ = I.saturation(S)
t1 = time.time()
print(f"  Done ({t1-t0:.2f}s)")

# Compute quotient dimension
print("Computing quotient dimension (= number of solutions)...")
try:
    Q = R.quotient(I_sat)
    dim = Q.vector_space_dimension()
    print(f"\n  QUOTIENT DIMENSION = {dim}")
    print()
    if dim == 24:
        print("  [CONFIRMED] Total solution count is 24")
    else:
        print(f"  [WARNING] Expected 24, got {dim}")
except Exception as e:
    print(f"  Error: {e}")
    # Alternative: compute Groebner basis degree
    print("  Trying alternative method...")
    gb = I_sat.groebner_basis()
    print(f"  Groebner basis has {len(gb)} elements")
