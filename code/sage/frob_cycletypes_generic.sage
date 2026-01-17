#!/usr/bin/env sage
"""
Frobenius Cycle Types via Eliminant with Generic 8D/11D Kinematics
===================================================================

This script extends frob_cycletypes_eliminant.sage to use generic kinematics
from 8D/11D lightcone generators, testing whether the Galois structure
persists when we move away from 4D special loci.

The key question: Is the D₂₄ Galois group and χ₈⁻ pattern a genuine feature
of the CHY cover, or an artifact of 4D kinematics hitting a special locus?

Usage:
    sage frob_cycletypes_generic.sage
    sage frob_cycletypes_generic.sage --kin_mode mandelstam --p_max 200
    sage frob_cycletypes_generic.sage --kin_mode 8d11d --dim 11 --p_max 100

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

# Load the new 8D/11D kinematics generators
load("generic_kinematics_8d11d.sage")

# Also load existing 4D for comparison
load("kinematics.sage")

print("=" * 70)
print("FROBENIUS CYCLE TYPES WITH GENERIC KINEMATICS")
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

kin_mode = _get_arg("--kin_mode", "mandelstam")  # 4d, 8d11d, mandelstam
dim = int(_get_arg("--dim", "11"))
n = int(_get_arg("--n", "7"))
p_min = int(_get_arg("--p_min", "5"))
p_max = int(_get_arg("--p_max", "200"))
seed = int(_get_arg("--seed", "0"))
outdir = _get_arg("--outdir", "../../logs")
gauge_str = _get_arg("--gauge", "0,1,-1")

gauge = tuple(int(x.strip()) for x in gauge_str.split(","))

print(f"Kinematics mode: {kin_mode}")
print(f"n = {n}")
if kin_mode == '8d11d':
    print(f"Dimension D = {dim}")
print(f"Prime range: [{p_min}, {p_max}]")
print(f"Seed: {seed}")
print(f"Gauge: {gauge}")
print()


# ============================================================================
# Kinematics generation
# ============================================================================

def generate_kin_for_eliminant(kin_mode, dim, n, seed):
    """
    Generate kinematics over QQ for later mod-p reduction.
    
    For 4d and 8d11d modes, we generate over QQ and reduce later.
    For mandelstam mode, we generate fresh for each prime.
    
    Returns:
        (mandelstams_dict, mode_is_per_prime)
    """
    if kin_mode == '4d':
        # Existing 4D generator over QQ
        momenta = make_4d_massless_point(n, seed)
        sij = momenta_to_mandelstams(momenta)
        return sij, False
    
    elif kin_mode == 'mandelstam':
        # Per-prime random generation
        return None, True
    
    elif kin_mode == '8d11d':
        # Generate over QQ is tricky for lightcone...
        # We'll use per-prime generation instead
        return None, True
    
    else:
        raise ValueError(f"Unknown kinematics mode: {kin_mode}")


def get_mandelstams_for_prime(kin_mode, dim, n, p, seed, base_mandelstams=None):
    """
    Get Mandelstams over GF(p) for eliminant computation.
    """
    F = GF(p)
    
    if base_mandelstams is not None:
        # Reduce QQ kinematics mod p
        sij_p = {}
        for key, val in base_mandelstams.items():
            sij_p[key] = F(val)
        return sij_p
    
    elif kin_mode == 'mandelstam':
        # Dimension-free generator
        _, sij_matrix = random_chy_mandelstams(n, p, use_p2=False, 
                                               avoid_zero=True, seed=seed + p)
        sij = {}
        for i in range(n):
            for j in range(i + 1, n):
                sij[(i + 1, j + 1)] = sij_matrix[i, j]
        return sij
    
    elif kin_mode == '8d11d':
        # 8D/11D generator
        try:
            _, momenta, sij_matrix = generate_null_kinematics_GFp(dim, n, p, 
                                                                   avoid_zero_sij=True,
                                                                   seed=seed + p)
            sij = {}
            for i in range(n):
                for j in range(i + 1, n):
                    sij[(i + 1, j + 1)] = sij_matrix[i, j]
            return sij
        except RuntimeError:
            return None
    
    else:
        raise ValueError(f"Unknown kinematics mode: {kin_mode}")


# ============================================================================
# CHY ideal and eliminant (adapted from frob_cycletypes_eliminant.sage)
# ============================================================================

def build_chy_ideal_over_fp(mandelstams, p, gauge):
    """
    Build CHY ideal over GF(p)[s7, s6, s5, s4] with lex order.
    """
    F = GF(p)
    
    # Lex order for elimination: s7 > s6 > s5 > s4
    R = PolynomialRing(F, names=["s7", "s6", "s5", "s4"], order="lex")
    s7, s6, s5, s4 = R.gens()
    
    # Gauge-fixed sigmas
    s1, s2, s3 = [F(g) for g in gauge]
    sig = {1: s1, 2: s2, 3: s3, 4: s4, 5: s5, 6: s6, 7: s7}
    
    def s_ab(a, b):
        i, j = min(a, b), max(a, b)
        return F(mandelstams[(i, j)])
    
    polys = []
    for a in [4, 5, 6, 7]:
        # denom_a = ∏_{b!=a} (σ_a - σ_b)
        denom_a = R(1)
        for b in range(1, 8):
            if b != a:
                denom_a *= (sig[a] - sig[b])
        
        # Cleared numerator
        num = R(0)
        for b in range(1, 8):
            if b != a:
                term = s_ab(a, b)
                num += term * (denom_a // (sig[a] - sig[b]))
        
        polys.append(num)
    
    I = Ideal(polys)
    return R, I


def eliminant_in_s4(R, I):
    """
    Extract univariate eliminant in s4 by eliminating s7, s6, s5.
    """
    s7, s6, s5, s4 = R.gens()
    
    J = I.elimination_ideal([s7, s6, s5])
    gens = J.gens()
    
    univs = [g for g in gens if set(g.variables()) <= set([s4]) and g != 0]
    
    if not univs:
        return None
    
    f = min(univs, key=lambda h: h.degree())
    
    try:
        f = f.monic()
    except Exception:
        pass
    
    return f


def is_squarefree_poly(f):
    if f.degree() <= 0:
        return False
    return gcd(f, f.derivative()) == 1


def cycle_type_from_factor_degrees(f):
    facs = f.factor()
    degs = []
    for g, e in facs:
        degs += [g.degree()] * e
    degs.sort(reverse=True)
    return degs


def classify_D24_rotation(degs):
    """Check if degrees match D24 rotation templates."""
    degs_sorted = sorted(degs)
    templates = [
        [1] * 24,
        [2] * 12,
        [3] * 8,
        [4] * 6,
        [6] * 4,
        [8] * 3,
        [12] * 2,
        [24] * 1,
    ]
    return degs_sorted in templates


# ============================================================================
# Main sweep
# ============================================================================

def main():
    os.makedirs(outdir, exist_ok=True)
    
    # Get base kinematics if using 4d mode
    base_mandelstams, per_prime = generate_kin_for_eliminant(kin_mode, dim, n, seed)
    
    results = {
        'timestamp': str(datetime.now()),
        'kin_mode': kin_mode,
        'dim': dim if kin_mode == '8d11d' else None,
        'n': n,
        'p_range': [p_min, p_max],
        'seed': seed,
        'gauge': list(gauge),
        'records': [],
        'bad_reduction': [],
        'summary': {
            'D24_template_matches': 0,
            'total_good': 0,
            'cycle_type_distribution': {}
        }
    }
    
    print("-" * 70)
    print(f"{'p':>5s} {'mod8':>5s} {'χ₄':>4s} {'χ₈⁻':>4s} {'deg':>4s} {'sq':>3s} {'cycle_type':>25s} {'D24':>4s} {'time':>7s}")
    print("-" * 70)
    
    total_time = 0
    
    for p in prime_range(p_min, p_max + 1):
        if p < 5:
            continue
        
        t0 = time.time()
        
        try:
            # Get Mandelstams
            mand = get_mandelstams_for_prime(kin_mode, dim, n, p, seed, base_mandelstams)
            if mand is None:
                results['bad_reduction'].append({'p': int(p), 'reason': 'kinematics_failed'})
                continue
            
            # Build ideal and extract eliminant
            R, I = build_chy_ideal_over_fp(mand, p, gauge)
            f = eliminant_in_s4(R, I)
            
            if f is None:
                results['bad_reduction'].append({'p': int(p), 'reason': 'no_eliminant'})
                continue
            
            deg = f.degree()
            if deg != 24:
                results['bad_reduction'].append({'p': int(p), 'reason': f'degree_{deg}'})
                continue
            
            sqfree = is_squarefree_poly(f)
            if not sqfree:
                results['bad_reduction'].append({'p': int(p), 'reason': 'not_squarefree'})
                continue
            
            # Get cycle type
            degs = cycle_type_from_factor_degrees(f)
            degs_int = [int(d) for d in degs]
            degs_tuple = tuple(degs_int)
            
            is_D24 = classify_D24_rotation(degs)
            
            elapsed = time.time() - t0
            total_time += elapsed
            
            c4 = chi4(p)
            c8m2 = chi8m2(p)
            c8p = chi8_2(p)
            
            rec = {
                'p': int(p),
                'mod4': int(p % 4),
                'mod8': int(p % 8),
                'chi4': int(c4),
                'chi8m2': int(c8m2),
                'chi8_2': int(c8p),
                'degree': int(deg),
                'squarefree': sqfree,
                'cycle_type': degs_int,
                'D24_template': is_D24,
                'time': round(elapsed, 3)
            }
            results['records'].append(rec)
            
            # Update summary
            results['summary']['total_good'] += 1
            if is_D24:
                results['summary']['D24_template_matches'] += 1
            
            cycle_key = str(degs_tuple)
            if cycle_key not in results['summary']['cycle_type_distribution']:
                results['summary']['cycle_type_distribution'][cycle_key] = {
                    'count': 0, 'primes': [], 'chi8_plus1': 0, 'chi8_minus1': 0
                }
            results['summary']['cycle_type_distribution'][cycle_key]['count'] += 1
            results['summary']['cycle_type_distribution'][cycle_key]['primes'].append(int(p))
            if c8m2 == 1:
                results['summary']['cycle_type_distribution'][cycle_key]['chi8_plus1'] += 1
            else:
                results['summary']['cycle_type_distribution'][cycle_key]['chi8_minus1'] += 1
            
            # Print
            d24_str = "✓" if is_D24 else "✗"
            sq_str = "Y" if sqfree else "N"
            cycle_str = str(degs_int)[:25]
            print(f"{p:5d} {p%8:5d} {c4:+4d} {c8m2:+4d} {deg:4d} {sq_str:>3s} {cycle_str:>25s} {d24_str:>4s} {elapsed:6.2f}s")
        
        except Exception as e:
            elapsed = time.time() - t0
            results['bad_reduction'].append({'p': int(p), 'reason': str(e)[:50]})
            print(f"{p:5d} {'---':>5s} {'---':>4s} {'---':>4s} {'---':>4s} {'---':>3s} {'FAILED':>25s} {'---':>4s} {elapsed:6.2f}s")
    
    # Summary
    print()
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    
    total = results['summary']['total_good']
    d24_matches = results['summary']['D24_template_matches']
    
    print(f"Total good primes: {total}")
    print(f"D₂₄ template matches: {d24_matches} ({100.0*d24_matches/total:.1f}%)" if total > 0 else "")
    print(f"Bad reduction: {len(results['bad_reduction'])}")
    print()
    
    print("Cycle type distribution:")
    for cycle_key, info in sorted(results['summary']['cycle_type_distribution'].items(), 
                                   key=lambda x: -x[1]['count']):
        n_plus = info['chi8_plus1']
        n_minus = info['chi8_minus1']
        print(f"  {cycle_key}: {info['count']} primes (χ₈⁻=+1: {n_plus}, χ₈⁻=-1: {n_minus})")
    
    # χ₈ correlation on cycle types
    print()
    print("-" * 70)
    print("χ₈⁻ CORRELATION ON CYCLE TYPES")
    print("-" * 70)
    
    # Check if certain cycle types appear more for χ₈ = +1 vs -1
    for cycle_key, info in sorted(results['summary']['cycle_type_distribution'].items(),
                                   key=lambda x: -x[1]['count']):
        if info['count'] >= 5:  # Only analyze if enough data
            total_ct = info['count']
            plus1_frac = info['chi8_plus1'] / total_ct
            minus1_frac = info['chi8_minus1'] / total_ct
            
            # Expected: roughly 50/50 if no correlation
            if abs(plus1_frac - 0.5) > 0.2:
                print(f"  {cycle_key}: χ₈⁻=+1 fraction = {100.0*plus1_frac:.0f}% (potential correlation)")
            else:
                print(f"  {cycle_key}: χ₈⁻=+1 fraction = {100.0*plus1_frac:.0f}% (no significant bias)")
    
    # Save results
    outfile = os.path.join(outdir, f"frob_generic_{kin_mode}_n{n}_seed{seed}.json")
    with open(outfile, "w") as f:
        json.dump(results, f, indent=2, cls=SageJSONEncoder)
    
    print()
    print("=" * 70)
    print(f"Results saved to: {outfile}")
    print(f"Total time: {total_time:.1f}s")
    print("=" * 70)


if __name__ == "__main__":
    main()
