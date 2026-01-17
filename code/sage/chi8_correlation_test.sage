#!/usr/bin/env sage
"""
χ₈ Correlation Test for CHY Scattering Equations
=================================================

This script tests whether the CHY solution count N_p correlates with χ₈⁻(p).

Key hypothesis being tested:
- If Frobenius splitting has a dyadic quadratic twist, the cycle structure
  should depend on χ₈⁻(p) = (-2/p).

Three kinematic generation modes:
1. '4d' - Standard 4D Pythagorean kinematics (existing code)
2. '8d11d' - Generic 8D/11D lightcone kinematics (new)
3. 'mandelstam' - Dimension-free CHY Mandelstams (new, fastest)

Usage:
    sage chi8_correlation_test.sage
    sage chi8_correlation_test.sage --mode mandelstam --p_max 200
    sage chi8_correlation_test.sage --mode 8d11d --dim 11 --p_max 100

Author: CHY arithmetic verification
Date: 2026-01-15
"""

from sage.all import *
import argparse
import json
import os
import sys
import time
from datetime import datetime


class SageJSONEncoder(json.JSONEncoder):
    """JSON encoder that handles Sage types."""
    def default(self, obj):
        # Handle Sage Integer
        if hasattr(obj, '__index__'):
            return int(obj)
        # Handle Sage Rational and RealDoubleElement
        try:
            return float(obj)
        except (TypeError, ValueError):
            pass
        # Fallback: convert to string
        try:
            return str(obj)
        except:
            pass
        return super().default(obj)

# Get the directory of this script for relative imports
_script_dir = os.path.dirname(os.path.abspath(__file__)) if '__file__' in dir() else os.getcwd()
os.chdir(_script_dir)

# Load our new 8D/11D generator
load("generic_kinematics_8d11d.sage")

# Also load existing 4D kinematics for comparison
load("kinematics.sage")

print("=" * 70)
print("χ₈⁻ CORRELATION TEST FOR CHY SCATTERING EQUATIONS")
print("=" * 70)
print(f"Timestamp: {datetime.now()}")
print()

# ============================================================================
# CLI argument parsing
# ============================================================================

def _get_arg(flag, default=None):
    if flag in sys.argv:
        i = sys.argv.index(flag)
        if i + 1 < len(sys.argv):
            return sys.argv[i + 1]
    return default

def _has_flag(flag):
    return flag in sys.argv

mode = _get_arg("--mode", "mandelstam")  # 4d, 8d11d, mandelstam
dim = int(_get_arg("--dim", "11"))
n = int(_get_arg("--n", "7"))
p_min = int(_get_arg("--p_min", "5"))
p_max = int(_get_arg("--p_max", "100"))
seed = int(_get_arg("--seed", "0"))
outdir = _get_arg("--outdir", "../../logs")
verbose = not _has_flag("--quiet")

gauge = (0, 1, -1)

print(f"Mode: {mode}")
print(f"n = {n}")
if mode == '8d11d':
    print(f"Dimension D = {dim}")
print(f"Prime range: [{p_min}, {p_max}]")
print(f"Seed: {seed}")
print(f"Gauge: {gauge}")
print()

# ============================================================================
# Kinematics generation dispatch
# ============================================================================

def generate_kinematics_for_prime(p, mode, dim, n, seed):
    """
    Generate kinematics over GF(p) for prime p.
    
    Returns:
        (sij_dict or sij_matrix, field)
    """
    if mode == '4d':
        # Existing 4D generator
        momenta = make_4d_massless_point(n, seed)
        sij_dict = momenta_to_mandelstams(momenta)
        return sij_dict, GF(p)
    
    elif mode == '8d11d':
        # New 8D/11D generator
        F, momenta, sij_matrix = generate_null_kinematics_GFp(dim, n, p, seed=seed)
        # Convert matrix to dict for consistency
        sij_dict = {}
        for i in range(n):
            for j in range(i + 1, n):
                sij_dict[(i + 1, j + 1)] = sij_matrix[i, j]
        return sij_dict, F
    
    elif mode == 'mandelstam':
        # Dimension-free generator
        F, sij_matrix = random_chy_mandelstams(n, p, use_p2=False, seed=seed)
        sij_dict = {}
        for i in range(n):
            for j in range(i + 1, n):
                sij_dict[(i + 1, j + 1)] = sij_matrix[i, j]
        return sij_dict, F
    
    else:
        raise ValueError(f"Unknown mode: {mode}")


def count_chy_solutions_from_dict(sij_dict, n, F, gauge):
    """
    Count CHY solutions from a dictionary of Mandelstams.
    """
    sigma1, sigma2, sigma3 = gauge
    
    var_names = [f'x{i}' for i in range(4, n + 1)]
    R = PolynomialRing(F, var_names)
    xs = R.gens()
    
    sigmas = {1: F(sigma1), 2: F(sigma2), 3: F(sigma3)}
    for idx, x in enumerate(xs):
        sigmas[idx + 4] = x
    
    def get_s(i, j):
        key = (min(i, j), max(i, j))
        val = sij_dict.get(key, sij_dict.get((i, j), 0))
        return F(val)
    
    # Build cleared CHY polynomials
    polys = []
    for a in range(4, n + 1):
        h_cleared = R(0)
        for b in range(1, n + 1):
            if b != a:
                term = get_s(a, b)
                for c in range(1, n + 1):
                    if c != a and c != b:
                        term *= (sigmas[a] - sigmas[c])
                h_cleared += term
        polys.append(h_cleared)
    
    I = R.ideal(polys)
    
    # Spurious factor
    S = R(1)
    for x in xs:
        S *= x * (x - F(sigma1)) * (x - F(sigma2)) * (x - F(sigma3))
    for i in range(len(xs)):
        for j in range(i + 1, len(xs)):
            S *= (xs[i] - xs[j])
    
    # Saturate and count
    I_sat, _ = I.saturation(S)
    
    try:
        sols = I_sat.variety()
        return len(sols)
    except Exception:
        return -1


# ============================================================================
# Main sweep
# ============================================================================

def main():
    os.makedirs(outdir, exist_ok=True)
    
    results = {
        'timestamp': str(datetime.now()),
        'mode': mode,
        'dim': dim if mode == '8d11d' else None,
        'n': n,
        'p_range': [p_min, p_max],
        'seed': seed,
        'gauge': list(gauge),
        'records': [],
        'summary': {
            'chi8_plus1': {'primes': [], 'counts': []},
            'chi8_minus1': {'primes': [], 'counts': []},
            'chi4_plus1': {'primes': [], 'counts': []},
            'chi4_minus1': {'primes': [], 'counts': []},
        }
    }
    
    expected_count = factorial(n - 3)  # (n-3)! = 24 for n=7
    
    print("-" * 70)
    print(f"Sweeping primes {p_min} to {p_max}...")
    print(f"Expected solution count: {expected_count}")
    print("-" * 70)
print(f"{'p':>5s} {'mod8':>5s} {'χ₄':>4s} {'χ₈⁻':>4s} {'N':>5s} {'status':>10s} {'time':>8s}")
    print("-" * 70)
    
    total_time = 0
    
    for p in prime_range(p_min, p_max + 1):
        if p < 5:
            continue
        
        t0 = time.time()
        
        try:
            # Generate kinematics
            sij, F = generate_kinematics_for_prime(p, mode, dim, n, seed + p)
            
            # Count solutions
            count = count_chy_solutions_from_dict(sij, n, F, gauge)
            
            elapsed = time.time() - t0
            total_time += elapsed
            
            # Characters
            c4 = chi4(p)
            c8m2 = chi8m2(p)
            c8p = chi8_2(p)
            
            # Determine status
            if count == expected_count:
                status = "OK"
            elif count == 0:
                status = "VANISH"
            elif count > 0:
                status = f"N={count}"
            else:
                status = "ERROR"
            
            rec = {
                'p': int(p),
                'mod4': int(p % 4),
                'mod8': int(p % 8),
                'chi4': int(c4),
                'chi8m2': int(c8m2),
                'chi8_2': int(c8p),
                'N': int(count) if count >= 0 else None,
                'time': round(elapsed, 3),
                'status': status
            }
            results['records'].append(rec)
            
            # Group by characters
            if c8m2 == 1:
                results['summary']['chi8_plus1']['primes'].append(int(p))
                results['summary']['chi8_plus1']['counts'].append(int(count))
            else:
                results['summary']['chi8_minus1']['primes'].append(int(p))
                results['summary']['chi8_minus1']['counts'].append(int(count))
            
            if c4 == 1:
                results['summary']['chi4_plus1']['primes'].append(int(p))
                results['summary']['chi4_plus1']['counts'].append(int(count))
            else:
                results['summary']['chi4_minus1']['primes'].append(int(p))
                results['summary']['chi4_minus1']['counts'].append(int(count))
            
            if verbose:
                print(f"{p:5d} {p%8:5d} {c4:+4d} {c8m2:+4d} {count:5d} {status:>10s} {elapsed:7.2f}s")
        
        except Exception as e:
            elapsed = time.time() - t0
            if verbose:
                print(f"{p:5d} {'---':>5s} {'---':>4s} {'---':>4s} {'---':>5s} {'FAILED':>10s} {elapsed:7.2f}s")
                print(f"       Error: {e}")
    
    # Summary statistics
    print()
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    
    for key in ['chi8_plus1', 'chi8_minus1', 'chi4_plus1', 'chi4_minus1']:
        counts = results['summary'][key]['counts']
        if counts:
            n_primes = len(counts)
            n_24 = sum(1 for c in counts if c == expected_count)
            n_0 = sum(1 for c in counts if c == 0)
            n_other = n_primes - n_24 - n_0
            
            results['summary'][key]['n_primes'] = int(n_primes)
            results['summary'][key]['n_expected'] = int(n_24)
            results['summary'][key]['n_vanish'] = int(n_0)
            results['summary'][key]['n_other'] = int(n_other)
            results['summary'][key]['mean'] = float(sum(counts)) / n_primes
            results['summary'][key]['distinct_values'] = [int(x) for x in sorted(set(counts))]
            
            label = key.replace('_', ' = ').replace('plus', '+').replace('minus', '-')
            print(f"{label}:")
            print(f"  Primes tested: {n_primes}")
            print(f"  N = {expected_count}: {n_24} ({100.0*n_24/n_primes:.1f}%)")
            print(f"  N = 0: {n_0} ({100.0*n_0/n_primes:.1f}%)")
            print(f"  Other: {n_other}")
            print(f"  Mean N: {float(sum(counts))/n_primes:.2f}")
            print(f"  Distinct values: {sorted(set(counts))}")
            print()
    
    # χ₈ correlation analysis
    print("-" * 70)
print("χ₈⁻ CORRELATION ANALYSIS")
    print("-" * 70)
    
    plus1_counts = results['summary']['chi8_plus1']['counts']
    minus1_counts = results['summary']['chi8_minus1']['counts']
    
    if plus1_counts and minus1_counts:
        # Check if there's a systematic difference
        plus1_mean = sum(plus1_counts) / len(plus1_counts)
        minus1_mean = sum(minus1_counts) / len(minus1_counts)
        
        plus1_24_rate = sum(1 for c in plus1_counts if c == expected_count) / len(plus1_counts)
        minus1_24_rate = sum(1 for c in minus1_counts if c == expected_count) / len(minus1_counts)
        
        plus1_0_rate = sum(1 for c in plus1_counts if c == 0) / len(plus1_counts)
        minus1_0_rate = sum(1 for c in minus1_counts if c == 0) / len(minus1_counts)
        
        print(f"χ₈⁻ = +1: mean N = {float(plus1_mean):.2f}, N={expected_count} rate = {100.0*plus1_24_rate:.1f}%, N=0 rate = {100.0*plus1_0_rate:.1f}%")
        print(f"χ₈⁻ = -1: mean N = {float(minus1_mean):.2f}, N={expected_count} rate = {100.0*minus1_24_rate:.1f}%, N=0 rate = {100.0*minus1_0_rate:.1f}%")
        
        # Verdict
        if abs(plus1_mean - minus1_mean) < 0.5 and abs(plus1_24_rate - minus1_24_rate) < 0.1:
            print("\n[RESULT] NO SIGNIFICANT χ₈⁻ CORRELATION DETECTED")
            print("         The solution count appears independent of χ₈⁻(p).")
            results['verdict'] = 'no_correlation'
        else:
            print(f"\n[RESULT] POTENTIAL χ₈⁻ CORRELATION DETECTED")
            print(f"         Mean difference: {plus1_mean - minus1_mean:.2f}")
            print(f"         Rate difference: {100.0*(plus1_24_rate - minus1_24_rate):.1f}%")
            results['verdict'] = 'correlation_detected'
    
    # Also check χ₄ (inert vs split)
    print()
    print("-" * 70)
    print("χ₄ CORRELATION ANALYSIS (inert vs split)")
    print("-" * 70)
    
    split_counts = results['summary']['chi4_plus1']['counts']
    inert_counts = results['summary']['chi4_minus1']['counts']
    
    if split_counts and inert_counts:
        split_mean = sum(split_counts) / len(split_counts)
        inert_mean = sum(inert_counts) / len(inert_counts)
        
        split_24_rate = sum(1 for c in split_counts if c == expected_count) / len(split_counts)
        inert_24_rate = sum(1 for c in inert_counts if c == expected_count) / len(inert_counts)
        
        split_0_rate = sum(1 for c in split_counts if c == 0) / len(split_counts)
        inert_0_rate = sum(1 for c in inert_counts if c == 0) / len(inert_counts)
        
        print(f"Split (p≡1 mod4): mean N = {float(split_mean):.2f}, N={expected_count} rate = {100.0*split_24_rate:.1f}%, N=0 rate = {100.0*split_0_rate:.1f}%")
        print(f"Inert (p≡3 mod4): mean N = {float(inert_mean):.2f}, N={expected_count} rate = {100.0*inert_24_rate:.1f}%, N=0 rate = {100.0*inert_0_rate:.1f}%")
        
        if inert_0_rate > 0.5:
            print("\n[RESULT] STRONG INERT VANISHING DETECTED")
            print(f"         {100.0*inert_0_rate:.0f}% of inert primes have N=0")
            results['inert_vanishing'] = True
        else:
            print("\n[RESULT] Inert vanishing not consistently observed")
            results['inert_vanishing'] = False
    
    # Save results
    outfile = os.path.join(outdir, f"chi8_test_{mode}_n{n}_seed{seed}.json")
    with open(outfile, "w") as f:
        json.dump(results, f, indent=2, cls=SageJSONEncoder)
    
    print()
    print("=" * 70)
    print(f"Results saved to: {outfile}")
    print(f"Total time: {total_time:.1f}s")
    print("=" * 70)


if __name__ == "__main__":
    main()
