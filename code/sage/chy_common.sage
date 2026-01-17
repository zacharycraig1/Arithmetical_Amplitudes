"""
Shared utilities for CHY scattering equations (tree and loop levels).
Includes support for rational kinematics.
"""

from sage.all import *

def _random_massless_vector_rational(scale=100):
    """Return (E, px, py, pz) in Q satisfying E^2 = px^2 + py^2 + pz^2."""
    while True:
        m = randint(1, scale)
        n = randint(1, scale)
        p = randint(1, scale)
        q = randint(1, scale)
        
        x = m^2 + n^2 - p^2 - q^2
        y = 2*(m*q + n*p)
        z = 2*(n*q - m*p)
        E = m^2 + n^2 + p^2 + q^2
        
        return (QQ(E)*choice([-1,1]), QQ(x)*choice([-1,1]), QQ(y)*choice([-1,1]), QQ(z)*choice([-1,1]))

def _massless_momenta_rational(n, scale=100):
    """Sample n massless momenta in Q^4 with total momentum zero using null completion."""
    while True:
        ks = [_random_massless_vector_rational(scale) for _ in range(n-2)]
        K = [sum(k[i] for k in ks) for i in range(4)]
        K_sq = K[0]^2 - K[1]^2 - K[2]^2 - K[3]^2
        
        if K_sq == 0: continue
        
        a = K[1]^2 - K[0]^2
        if a == 0: continue
        
        for _ in range(100):
            y = QQ(randint(-scale, scale))
            z = QQ(randint(-scale, scale))
            
            # (K1^2 - K0^2)*x^2 + 2*K1*C*x + (C^2 - K0^2*(y^2+z^2)) = 0
            # C = y*K2 + z*K3 - K_sq/2
            C = y*K[2] + z*K[3] - K_sq/2
            
            b = 2*K[1]*C
            c_val = C^2 - K[0]^2*(y^2+z^2)
            
            disc = b^2 - 4*a*c_val
            if disc >= 0 and disc.is_square():
                x = (-b + disc.sqrt()) / (2*a)
                E = (x*K[1] + y*K[2] + z*K[3] - K_sq/2) / K[0]
                
                kn_1 = (E, x, y, z)
                kn = tuple(-K[i] - kn_1[i] for i in range(4))
                
                res = ks + [kn_1, kn]
                return res

def _get_rational_kinematics(n, seed=42):
    """Generate or retrieve a fixed set of rational kinematics."""
    set_random_seed(seed)
    return _massless_momenta_rational(n), None

def _reduce_kinematics_mod_p(momenta_q, p):
    """Reduce rational kinematics modulo p."""
    F = GF(p)
    momenta_p = []
    for m in momenta_q:
        momenta_p.append(tuple(F(x) for x in m))
        
    invariants_p = {}
    n = len(momenta_p)
    for i in range(n):
        Ei, pxi, pyi, pzi = momenta_p[i]
        for j in range(i + 1, n):
            Ej, pxj, pyj, pzj = momenta_p[j]
            dot = Ei * Ej - pxi * pxj - pyi * pyj - pzi * pzj
            invariants_p[(i + 1, j + 1)] = 2 * dot
            
    return invariants_p, momenta_p

def _get_invariant(invariants, i, j):
    if i < j:
        return invariants[(i, j)]
    return invariants[(j, i)]

def _gauge_fixed_sigma(n, ring):
    """Return a mapping of punctures to sigma coordinates with gauge fixing."""
    fixed = {1: ring(0), 2: ring(1), 3: ring(-1)}
    symbolic = {
        i: ring.gen(i - 4)
        for i in range(4, n + 1)
    }
    fixed.update(symbolic)
    return fixed, list(symbolic.values())

def _collision_discriminant(sigmas, labels):
    """Compute Δ = ∏_{A<B} (σ_A - σ_B) for the provided labels."""
    if not labels:
        first = next(iter(sigmas.values()))
        return first.parent()(1)
    delta = sigmas[labels[0]].parent()(1)
    for idx, a in enumerate(labels):
        for b in labels[idx + 1:]:
            delta *= sigmas[a] - sigmas[b]
    return delta
