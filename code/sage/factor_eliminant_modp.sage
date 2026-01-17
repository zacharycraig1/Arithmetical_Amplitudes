#!/usr/bin/env sage
"""
INSANE GUT CHECK: Factor the eliminant polynomial mod p
See if it splits into 12+12 or reveals the factorization structure
"""

import sys
import time

def generate_4D_null_momenta(n, seed=0):
    """Generate n null momenta in 4D that sum to zero"""
    set_random_seed(seed)
    
    momenta = []
    for i in range(n-1):
        # Random spinors
        la = vector([QQ.random_element(num_bound=5), QQ.random_element(num_bound=5)])
        if la == vector([0,0]):
            la = vector([1, 1])
        mu = vector([QQ.random_element(num_bound=5), QQ.random_element(num_bound=5)])
        if mu == vector([0,0]):
            mu = vector([1, 1])
        
        # p^{aa'} = lambda^a * mu^{a'} (null by construction)
        p = [la[0]*mu[0], la[0]*mu[1], la[1]*mu[0], la[1]*mu[1]]
        momenta.append(vector(QQ, p))
    
    # Last momentum from conservation
    pn = -sum(momenta)
    momenta.append(pn)
    
    return momenta

def momenta_to_mandelstams(momenta):
    """Convert momenta to Mandelstam invariants s_ij = 2 p_i . p_j"""
    n = len(momenta)
    eta = diagonal_matrix([1, -1, -1, 1])  # Minkowski metric in spinor form
    
    s = {}
    for i in range(n):
        for j in range(i+1, n):
            s[i,j] = 2 * momenta[i] * eta * momenta[j]
            s[j,i] = s[i,j]
    return s

def build_CHY_ideal(n, s, R, z, gauge='standard'):
    """Build CHY scattering equations as polynomial ideal"""
    polys = []
    
    for i in range(1, n-1):  # Scattering equations for z_2, ..., z_{n-1}
        if gauge == 'standard' and (i == 0 or i == n-2 or i == n-1):
            continue
            
        fi = R(0)
        for j in range(n):
            if j != i:
                sij = s.get((i,j), s.get((j,i), 0))
                if j == 0:
                    zj = R(0)
                elif j == n-1:
                    continue  # z_n = infinity, handled separately
                elif j == n-2:
                    zj = R(1)
                else:
                    zj = z[j-1]  # z[0] = z_2, z[1] = z_3, etc.
                
                zi = z[i-1] if i > 0 else R(0)
                if j < n-1:
                    fi += sij / (zi - zj)
        
        # Clear denominators
        denom = R(1)
        for j in range(n-1):
            if j != i:
                if j == 0:
                    zj = R(0)
                elif j == n-2:
                    zj = R(1)
                else:
                    zj = z[j-1]
                zi = z[i-1]
                denom *= (zi - zj)
        
        poly = R(0)
        for j in range(n):
            if j != i and j < n-1:
                sij = s.get((i,j), s.get((j,i), 0))
                if j == 0:
                    zj = R(0)
                elif j == n-2:
                    zj = R(1)
                else:
                    zj = z[j-1]
                
                zi = z[i-1]
                term = sij
                for k in range(n-1):
                    if k != i and k != j:
                        if k == 0:
                            zk = R(0)
                        elif k == n-2:
                            zk = R(1)
                        else:
                            zk = z[k-1]
                        term *= (zi - zk)
                poly += term
        
        if poly != 0:
            polys.append(poly)
    
    return polys

def compute_eliminant_modp(n, s, p, var_idx=0):
    """
    Compute the eliminant polynomial for z_{var_idx+2} modulo p
    by elimination from the CHY ideal
    """
    F = GF(p)
    
    # Map mandelstams to F
    s_modp = {}
    for key, val in s.items():
        s_modp[key] = F(val)
    
    # Create polynomial ring over F
    n_vars = n - 3  # z_2, z_3, ..., z_{n-2}
    var_names = ['z%d' % (i+2) for i in range(n_vars)]
    R = PolynomialRing(F, var_names, order='lex')
    z = R.gens()
    
    # Build ideal
    polys = []
    
    # Scattering equations
    for a in range(1, n-2):  # a = 1, ..., n-3 (for z_2, ..., z_{n-2})
        eq = R(0)
        za = z[a-1]
        
        for b in range(n):
            if b == a:
                continue
            sab = s_modp.get((a,b), s_modp.get((b,a), F(0)))
            if sab == 0:
                continue
                
            if b == 0:
                zb = F(0)
            elif b == n-1:
                # Handle z_n = infinity: term is s_{an}/z_a
                eq += sab / za
                continue
            elif b == n-2:
                zb = F(1)
            else:
                zb = z[b-1]
            
            eq += sab / (za - zb)
        
        # Clear denominators
        factors = [za]  # from z_n = infinity term
        for b in range(n-1):
            if b == a:
                continue
            if b == 0:
                zb = F(0)
            elif b == n-2:
                zb = F(1)
            else:
                zb = z[b-1]
            factors.append(za - zb)
        
        common = prod(factors)
        poly = R(0)
        
        for b in range(n):
            if b == a:
                continue
            sab = s_modp.get((a,b), s_modp.get((b,a), F(0)))
            if sab == 0:
                continue
            
            if b == n-1:
                term = sab * prod(f for f in factors if f != za)
            else:
                if b == 0:
                    zb = F(0)
                elif b == n-2:
                    zb = F(1)
                else:
                    zb = z[b-1]
                term = sab * prod(f for f in factors if f != (za - zb))
            
            poly += term
        
        polys.append(poly)
    
    # Add non-collision constraints if needed
    
    print(f"Built {len(polys)} polynomials in {R}")
    for i, p in enumerate(polys):
        print(f"  f_{i+1}: degree {p.degree()}, {len(p.monomials())} terms")
    
    # Compute Groebner basis
    print(f"\nComputing Groebner basis over GF({p})...")
    t0 = time.time()
    I = R.ideal(polys)
    
    try:
        G = I.groebner_basis()
        t1 = time.time()
        print(f"Groebner basis computed in {t1-t0:.2f}s, {len(G)} elements")
        
        # Look for univariate polynomials
        target_var = z[var_idx]
        
        print(f"\nLooking for eliminant in {target_var}...")
        for g in G:
            vars_in_g = [v for v in z if g.degree(v) > 0]
            if len(vars_in_g) == 1 and vars_in_g[0] == target_var:
                print(f"\n*** FOUND ELIMINANT ***")
                print(f"Degree: {g.degree()}")
                print(f"Polynomial: {g}")
                
                # Factor it!
                print(f"\n*** FACTORING ***")
                fact = g.factor()
                print(f"Factorization: {fact}")
                
                degrees = []
                for f, e in fact:
                    d = f.degree()
                    degrees.extend([d] * e)
                print(f"Factor degrees: {sorted(degrees, reverse=True)}")
                print(f"Sum of degrees: {sum(degrees)} (should be 24)")
                
                return g, fact
        
        print("No univariate eliminant found in Groebner basis")
        print("Basis elements:")
        for g in G[:5]:
            print(f"  {g}")
        
        return None, None
        
    except Exception as e:
        print(f"Groebner basis failed: {e}")
        return None, None

def main():
    n = 7
    seed = int(sys.argv[1]) if len(sys.argv) > 1 else 42
    p = int(sys.argv[2]) if len(sys.argv) > 2 else 97
    
    print(f"=" * 60)
    print(f"ELIMINANT FACTORIZATION TEST")
    print(f"n={n}, seed={seed}, p={p}")
    print(f"=" * 60)
    
    # Generate kinematics
    print(f"\nGenerating kinematics with seed {seed}...")
    momenta = generate_4D_null_momenta(n, seed)
    s = momenta_to_mandelstams(momenta)
    
    # Check for zero mandelstams
    zero_count = sum(1 for v in s.values() if v == 0)
    print(f"Mandelstams: {len(s)} values, {zero_count} zeros")
    
    # Compute eliminant
    elim, fact = compute_eliminant_modp(n, s, p, var_idx=0)
    
    if elim is not None:
        print(f"\n" + "=" * 60)
        print(f"SUCCESS! Eliminant found and factored.")
        print(f"=" * 60)

if __name__ == "__main__":
    main()
