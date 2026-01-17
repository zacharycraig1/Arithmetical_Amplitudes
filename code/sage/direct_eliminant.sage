#!/usr/bin/env sage
"""
Direct eliminant computation via resultant cascade
"""

import sys
import time

def main():
    n = 7
    seed = int(sys.argv[1]) if len(sys.argv) > 1 else 42
    p = int(sys.argv[2]) if len(sys.argv) > 2 else 97
    
    print(f"=" * 60)
    print(f"DIRECT ELIMINANT VIA VARIETY")
    print(f"n={n}, seed={seed}, p={p}")
    print(f"=" * 60)
    
    # Use the existing orbit_size_detector approach which works
    # but extract z2 values instead of counting
    
    set_random_seed(seed)
    F = GF(p)
    
    # Generate random 4D null momenta
    momenta = []
    for i in range(n-1):
        la = vector([QQ.random_element(num_bound=5), QQ.random_element(num_bound=5)])
        if la == vector([0,0]):
            la = vector([1, 1])
        mu = vector([QQ.random_element(num_bound=5), QQ.random_element(num_bound=5)])
        if mu == vector([0,0]):
            mu = vector([1, 1])
        p_vec = [la[0]*mu[0], la[0]*mu[1], la[1]*mu[0], la[1]*mu[1]]
        momenta.append(vector(QQ, p_vec))
    pn = -sum(momenta)
    momenta.append(pn)
    
    # Mandelstams
    eta = diagonal_matrix([1, -1, -1, 1])
    s = {}
    for i in range(n):
        for j in range(i+1, n):
            s[i,j] = 2 * momenta[i] * eta * momenta[j]
            s[j,i] = s[i,j]
    
    print(f"Generated {len(s)} mandelstams")
    
    # Build ideal over F
    n_vars = n - 3
    var_names = ['z%d' % (i+2) for i in range(n_vars)]
    R = PolynomialRing(F, var_names)
    z = R.gens()
    
    polys = []
    for a in range(1, n-2):
        za = z[a-1]
        
        # Build cleared-denominator form of scattering equation
        factors = [za]  # from 1/z_a term when b=n
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
        
        poly = R(0)
        for b in range(n):
            if b == a:
                continue
            sab = F(s.get((a,b), s.get((b,a), 0)))
            if sab == 0:
                continue
            
            if b == n-1:  # b = 7, z_7 = infinity
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
        
        if poly != 0:
            polys.append(poly)
    
    print(f"Built {len(polys)} scattering equations")
    
    # Try to find solutions via variety
    print(f"\nComputing variety over F_{p}...")
    t0 = time.time()
    I = R.ideal(polys)
    
    # Saturate to remove non-collision
    non_coll = R(1)
    for i in range(n_vars):
        zi = z[i]
        # z_i != 0
        non_coll *= zi
        # z_i != 1
        non_coll *= (zi - 1)
        for j in range(i+1, n_vars):
            zj = z[j]
            non_coll *= (zi - zj)
    
    print(f"Saturating by non-collision constraints...")
    I_sat = I.saturation(non_coll)[0]
    
    try:
        V = I_sat.variety()
        t1 = time.time()
        print(f"Found {len(V)} solutions in {t1-t0:.2f}s")
        
        if len(V) > 0:
            # Extract z2 values
            z2_vals = [sol[z[0]] for sol in V]
            print(f"\nz2 values: {z2_vals}")
            
            # Build eliminant as product of (x - z2_i)
            Rx = PolynomialRing(F, 'x')
            x = Rx.gen()
            elim = prod(x - v for v in z2_vals)
            
            print(f"\nEliminant (in z2):")
            print(f"  Degree: {elim.degree()}")
            print(f"  Polynomial: {elim}")
            
            # Factor it
            print(f"\n*** FACTORING ELIMINANT ***")
            fact = elim.factor()
            print(f"Factorization: {fact}")
            
            degrees = []
            for f, e in fact:
                d = f.degree()
                degrees.extend([d] * e)
            degrees.sort(reverse=True)
            print(f"Factor degrees: {degrees}")
            print(f"Partition: {' + '.join(map(str, degrees))} = {sum(degrees)}")
            
            # Look for patterns
            if degrees.count(12) == 2 and sum(degrees) == 24:
                print(f"\n*** 12 + 12 STRUCTURE CONFIRMED! ***")
            elif max(degrees) == 24:
                print(f"\n*** IRREDUCIBLE (degree 24) ***")
            else:
                print(f"\n*** INTERESTING FACTORIZATION PATTERN ***")
                
    except Exception as e:
        print(f"Variety computation failed: {e}")

if __name__ == "__main__":
    main()
