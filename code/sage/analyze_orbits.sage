#!/usr/bin/env sage
"""
Analyze orbit decomposition from N_{p^k} data
"""

# Data from p=97, seed=42
data = {
    1: 3, 2: 7, 3: 6, 4: 15, 5: 8, 6: 10,
    7: 3, 8: 15, 9: 6, 10: 12, 11: 3, 12: 18,
    13: 3, 14: 7, 15: 11, 16: 15, 17: 3, 18: 10,
    19: 3, 20: 20, 21: 6, 22: 7, 23: 3, 24: 18
}

print("Analyzing Frobenius orbit structure from N_{p^k} data")
print("=" * 60)

# Count orbits by size
# N_{p^k} = sum over orbits O of 1_{d|k} where d = |O|
# So if we know N_{p^k} for all k, we can deduce orbit sizes.

# For orbit of size d, it contributes to N_{p^k} iff d | k

# Method: use Mobius inversion or direct counting
# Number of orbits of size exactly d = (N_{p^d} - N_{p^d-1}) / d 
# (approximately, need inclusion-exclusion)

# Actually, let's use: if orbit has size d, then d | k => contributes 1 per solution
# So N_{p^k} = sum_{d | k} (number of solutions in orbits of size d)

# Let a_d = number of solutions in orbits of size d
# Then N_{p^k} = sum_{d | k} a_d

# We can solve this via Mobius inversion:
# a_d = sum_{d | k, k | max_k} mu(k/d) * N_{p^k}

max_k = 24

def divisors(n):
    return [d for d in range(1, n+1) if n % d == 0]

def mobius(n):
    if n == 1:
        return 1
    factors = list(factor(n))
    if any(e > 1 for p, e in factors):
        return 0
    return (-1)**len(factors)

# For each divisor d of max_k, compute a_d
print("Computing orbit sizes using Mobius inversion...")
print()

orbit_counts = {}
for d in divisors(max_k):
    # a_d = sum_{d | k} mu(k/d) * N_{p^k}, where k ranges over multiples of d up to max_k
    a_d = 0
    for k in range(d, max_k + 1, d):
        m = mobius(k // d)
        if m != 0:
            a_d += m * data.get(k, 0)
    orbit_counts[d] = a_d
    if a_d != 0:
        print(f"  d = {d}: {a_d} solutions in orbits of size {d}")

print()
print("=" * 60)
print("ORBIT DECOMPOSITION:")
print("=" * 60)

total = 0
for d in sorted(orbit_counts.keys()):
    if orbit_counts[d] > 0:
        n_orbits = orbit_counts[d] // d
        print(f"  {n_orbits} orbits of size {d} (contributing {orbit_counts[d]} solutions)")
        total += orbit_counts[d]

print()
print(f"TOTAL SOLUTIONS: {total}")
print()

if total == 24:
    print("*** EXACTLY 24 SOLUTIONS - MATCHES CHY EXPECTATION ***")
elif total < 24:
    print(f"*** ONLY {total} SOLUTIONS VISIBLE! ***")
    print(f"*** {24 - total} SOLUTIONS MISSING (bad reduction or larger orbits) ***")
else:
    print(f"*** MORE THAN 24 SOLUTIONS - CHECK CALCULATION ***")

# Verify
print()
print("Verification:")
for k in [1, 2, 4, 8, 12, 24]:
    predicted = sum(orbit_counts[d] for d in divisors(k) if d in orbit_counts)
    actual = data.get(k, 0)
    match = "OK" if predicted == actual else f"MISMATCH (expected {predicted})"
    print(f"  N_{{p^{k}}} = {actual}  {match}")
