#!/usr/bin/env sage
"""
Fast eliminant factorization - skip saturation, filter solutions manually
"""

import sys
import time

def main():
    n = 7
    seed = int(sys.argv[1]) if len(sys.argv) > 1 else 42
    p = int(sys.argv[2]) if len(sys.argv) > 2 else 17  # Smaller prime
    
    print(f"=" * 60)
    print(f"FAST ELIMINANT TEST")
    print(f"n={n}, seed={seed}, p={p}")
    print(f"=" * 60)
    
    set_random_seed(seed)
    F = GF(p)
    
    # Generate momenta
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
    
    eta = diagonal_matrix([1, -1, -1, 1])
    s = {}
    for i in range(n):
        for j in range(i+1, n):
            s[i,j] = 2 * momenta[i] * eta * momenta[j]
            s[j,i] = s[i,j]
    
    # Build ideal
    n_vars = n - 3
    var_names = ['z%d' % (i+2) for i in range(n_vars)]
    R = PolynomialRing(F, var_names)
    z = R.gens()
    
    polys = []
    for a in range(1, n-2):
        za = z[a-1]
        factors = [za]
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
        
        if poly != 0:
            polys.append(poly)
    
    print(f"Built {len(polys)} equations")
    
    # Just find variety without saturation
    print(f"\nComputing variety...")
    t0 = time.time()
    I = R.ideal(polys)
    
    try:
        V = I.variety()
        t1 = time.time()
        print(f"Raw variety: {len(V)} solutions in {t1-t0:.2f}s")
        
        # Filter valid solutions (no collisions)
        valid = []
        for sol in V:
            vals = [sol[zi] for zi in z]
            # Check distinct and not 0 or 1
            points = [F(0)] + vals + [F(1)]  # z1=0, z2,...,z5, z6=1
            if len(set(points)) == len(points):  # all distinct
                valid.append(sol)
        
        print(f"Valid solutions (no collisions): {len(valid)}")
        
        if len(valid) > 0:
            z2_vals = [sol[z[0]] for sol in valid]
            print(f"z2 values: {z2_vals}")
            
            # Build eliminant
            Rx = PolynomialRing(F, 'x')
            x = Rx.gen()
            elim = prod(x - v for v in z2_vals)
            
            print(f"\nEliminant degree: {elim.degree()}")
            
            fact = elim.factor()
            print(f"\n*** FACTORIZATION ***")
            print(f"{fact}")
            
            degrees = []
            for f, e in fact:
                degrees.extend([f.degree()] * e)
            degrees.sort(reverse=True)
            print(f"\nFactor degrees: {degrees}")
            print(f"Sum: {sum(degrees)}")
            
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()
