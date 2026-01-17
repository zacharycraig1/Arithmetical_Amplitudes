#!/usr/bin/env sage
"""
Generic kinematics with momentum conservation via projection
"""

import sys
import time

print("=" * 60)
print("GENERIC KINEMATICS v3 (Projection Method)")
print("=" * 60)

def generate_generic_mandelstams(n, seed):
    """
    Generate s_ij satisfying momentum conservation using projection.
    
    Constraint: sum_{j!=i} s_ij = 0 for all i
    
    Method: Generate random s_ij, then subtract the constraint-violating part.
    """
    set_random_seed(seed)
    
    # Start with random values
    s = {}
    for i in range(1, n+1):
        for j in range(i+1, n+1):
            s[(i,j)] = QQ(randint(-50, 50))
            if s[(i,j)] == 0:
                s[(i,j)] = QQ(1)
    
    # Compute constraint violations
    violations = []
    for i in range(1, n+1):
        v = sum(s.get((min(i,j), max(i,j)), 0) for j in range(1, n+1) if j != i)
        violations.append(v)
    
    # Project: subtract (v_i + v_j)/(2n-2) from each s_ij
    # This ensures sum_j s'_ij = sum_j s_ij - sum_j (v_i + v_j)/(2n-2)
    #                         = v_i - (n-1)*v_i/(2n-2) - sum_j v_j/(2n-2)
    # Hmm, this gets complicated. Let me use a simpler method.
    
    # Simpler: Iteratively adjust until constraints are satisfied
    for iteration in range(100):
        max_violation = 0
        for i in range(1, n+1):
            v_i = sum(s.get((min(i,j), max(i,j)), 0) for j in range(1, n+1) if j != i)
            if abs(v_i) > abs(max_violation):
                max_violation = v_i
            
            if v_i != 0:
                # Distribute -v_i evenly among all s_ij for j != i
                adjustment = -v_i / (n - 1)
                for j in range(1, n+1):
                    if j != i:
                        key = (min(i,j), max(i,j))
                        s[key] = s[key] + adjustment
        
        if max_violation == 0:
            break
    
    # Verify
    ok = True
    for i in range(1, n+1):
        v_i = sum(s.get((min(i,j), max(i,j)), 0) for j in range(1, n+1) if j != i)
        if v_i != 0:
            print(f"  Residual at particle {i}: {v_i}")
            ok = False
    
    if ok:
        print("  Momentum conservation satisfied!")
    
    return s


def count_solutions(s, n, F, gauge):
    sigma1, sigma2, sigma3 = gauge
    
    var_names = ['x%d' % i for i in range(4, n+1)]
    R = PolynomialRing(F, var_names)
    xs = R.gens()
    
    sigmas = {1: F(sigma1), 2: F(sigma2), 3: F(sigma3)}
    for i, x in enumerate(xs):
        sigmas[i+4] = x
    
    def get_s(i, j):
        key = (min(i,j), max(i,j))
        return F(s[key])
    
    polys = []
    for a in range(4, n+1):
        h_cleared = R(0)
        for b in range(1, n+1):
            if b != a:
                term = get_s(a, b)
                for c in range(1, n+1):
                    if c != a and c != b:
                        term *= (sigmas[a] - sigmas[c])
                h_cleared += term
        polys.append(h_cleared)
    
    I = R.ideal(polys)
    
    S = R(1)
    for x in xs:
        S *= x * (x - F(sigma1)) * (x - F(sigma2)) * (x - F(sigma3))
    for i in range(len(xs)):
        for j in range(i+1, len(xs)):
            S *= (xs[i] - xs[j])
    
    I_sat, _ = I.saturation(S)
    sols = I_sat.variety()
    return len(sols)


n = 7
gauge = (0, 1, 2)
seed = int(sys.argv[1]) if len(sys.argv) > 1 else 0
p = int(sys.argv[2]) if len(sys.argv) > 2 else 101

print(f"n = {n}, seed = {seed}, p = {p}")
print()

s = generate_generic_mandelstams(n, seed)

print(f"\nCounting solutions over F_{p}...")
t0 = time.time()
count = count_solutions(s, n, GF(p), gauge)
t1 = time.time()

print(f"\n*** {count} SOLUTIONS over F_{p} ({t1-t0:.1f}s) ***")

if count == 24:
    print("\n***** EXACTLY 24! *****")
elif count > 20:
    print(f"\n*** CLOSE TO 24! ***")
