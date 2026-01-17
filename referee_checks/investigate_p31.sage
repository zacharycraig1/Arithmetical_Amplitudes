"""
INVESTIGATION: Why does p=31 fail for some seeds?
==================================================

This script performs the 5 fast checks to diagnose p=31 failures:
1. Denominator vanishing mod 31
2. s_ij = 0 mod 31 (soft/collinear degeneracy)
3. Collision saturation boundary solutions
4. Jacobian degeneracy (singular solution)
5. Discriminant/branch prime detection

Run: sage investigate_p31.sage
"""

import os
import sys
import json
from datetime import datetime

# Setup paths
_BASE_DIR = os.path.dirname(os.path.abspath(__file__))
_KINEMATICS_PATH = os.path.normpath(os.path.join(_BASE_DIR, "..", "code", "sage", "kinematics.sage"))
load(_KINEMATICS_PATH)

print("=" * 70)
print("INVESTIGATION: p=31 Failures")
print("=" * 70)
print(f"Timestamp: {datetime.now()}")
print()

# Configuration
P = 31
GAUGE = (0, 1, -1)
N = 7

# Load D1 seeds
def load_D1_seeds(path=None):
    if path is None:
        path = os.path.normpath(os.path.join(_BASE_DIR, "..", "data", "D1_seeds.json"))
    with open(path, "r") as f:
        obj = json.load(f)
    seeds = []
    for s in obj.get("seeds", []):
        if s.get("status", "valid") == "valid":
            seeds.append((s["id"], s["seed_value"]))
    return seeds

# CHY equation builder (same as main verification)
def build_cleared_polynomials(F, mandelstams, gauge):
    sigma1, sigma2, sigma3 = gauge
    R = PolynomialRing(F, 'x4,x5,x6,x7')
    x4, x5, x6, x7 = R.gens()
    sig = {1: F(sigma1), 2: F(sigma2), 3: F(sigma3), 4: x4, 5: x5, 6: x6, 7: x7}

    def s_ij(i, j):
        a, b = min(i, j), max(i, j)
        return F(mandelstams[(a, b)])

    polys = []
    for a in [4, 5, 6, 7]:
        expr = 0
        for b in [1, 2, 3, 4, 5, 6, 7]:
            if b == a:
                continue
            term = s_ij(a, b)
            for c in [1, 2, 3, 4, 5, 6, 7]:
                if c == a or c == b:
                    continue
                term *= (sig[a] - sig[c])
            expr += term
        polys.append(R(expr))
    return R, polys


def saturate_and_get_solutions(R, polys, gauge):
    """Return saturated ideal and explicit solutions."""
    sigma1, sigma2, sigma3 = gauge
    x4, x5, x6, x7 = R.gens()
    I = R.ideal(polys)

    S = 1
    for x in [x4, x5, x6, x7]:
        for g in [R.base_ring()(sigma1), R.base_ring()(sigma2), R.base_ring()(sigma3)]:
            S *= (x - g)
    S *= (x4 - x5) * (x4 - x6) * (x4 - x7) * (x5 - x6) * (x5 - x7) * (x6 - x7)

    I_sat, _ = I.saturation(S)
    
    try:
        sols = I_sat.variety()
    except:
        sols = []
    
    return I_sat, sols


def compute_chy_jacobian(mandelstams, sigma_vals, p):
    """
    Compute the CHY Jacobian determinant at a given solution.
    
    sigma_vals: dict {1: s1, 2: s2, ..., 7: s7}
    
    The Jacobian is d(h_4, h_5, h_6, h_7) / d(s_4, s_5, s_6, s_7)
    """
    F = GF(p)
    
    def s_ij(i, j):
        a, b = min(i, j), max(i, j)
        return F(mandelstams[(a, b)])
    
    # Partial derivatives: d h_a / d sigma_b where a, b in {4,5,6,7}
    # h_a = sum_{b != a} s_{ab} / (sigma_a - sigma_b)
    # d h_a / d sigma_a = - sum_{b != a} s_{ab} / (sigma_a - sigma_b)^2
    # d h_a / d sigma_b (b != a) = s_{ab} / (sigma_a - sigma_b)^2
    
    J = matrix(F, 4, 4)
    
    for row, a in enumerate([4, 5, 6, 7]):
        for col, b in enumerate([4, 5, 6, 7]):
            if a == b:
                # Diagonal: - sum_{c != a} s_{ac} / (sigma_a - sigma_c)^2
                total = F(0)
                for c in range(1, 8):
                    if c == a:
                        continue
                    denom = sigma_vals[a] - sigma_vals[c]
                    if denom == 0:
                        return None  # Collision
                    total -= s_ij(a, c) / (denom^2)
                J[row, col] = total
            else:
                # Off-diagonal: s_{ab} / (sigma_a - sigma_b)^2
                denom = sigma_vals[a] - sigma_vals[b]
                if denom == 0:
                    return None  # Collision
                J[row, col] = s_ij(a, b) / (denom^2)
    
    return J.det()


def investigate_seed(seed_id, seed_value, p):
    """Full investigation for a single seed at prime p."""
    print(f"\n{'='*70}")
    print(f"SEED {seed_id} (seed_value={seed_value}) at p={p}")
    print("="*70)
    
    # Generate kinematics
    momenta = make_4d_massless_point(N, seed_value)
    mandelstams = momenta_to_mandelstams(momenta)
    eps = levi_civita(momenta[0], momenta[1], momenta[2], momenta[3])
    
    F = GF(p)
    sigma1, sigma2, sigma3 = GAUGE
    
    findings = {
        "seed_id": seed_id,
        "seed_value": seed_value,
        "p": p,
        "issues": []
    }
    
    # =========================================================================
    # CHECK 1: Denominator vanishing mod p
    # =========================================================================
    print("\n[CHECK 1] Denominator vanishing mod p")
    print("-" * 40)
    
    denom_issues = []
    for key, val in mandelstams.items():
        if hasattr(val, 'denominator'):
            d = val.denominator()
            if d % p == 0:
                denom_issues.append((key, val, d))
                print(f"  s_{key} = {val} has denominator {d} ≡ 0 (mod {p})")
    
    # Also check epsilon
    if hasattr(eps, 'denominator') and eps.denominator() % p == 0:
        denom_issues.append(("epsilon", eps, eps.denominator()))
        print(f"  epsilon = {eps} has denominator ≡ 0 (mod {p})")
    
    if not denom_issues:
        print(f"  ✓ No denominators vanish mod {p}")
    else:
        findings["issues"].append({"type": "denominator", "details": str(denom_issues)})
    
    # =========================================================================
    # CHECK 2: s_ij = 0 mod p (soft/collinear degeneracy)
    # =========================================================================
    print("\n[CHECK 2] s_ij ≡ 0 (mod p)")
    print("-" * 40)
    
    # Skip if denominators vanish - this is bad reduction
    if denom_issues:
        print(f"\n  *** SKIPPING: Bad reduction (denominator vanishes mod {p}) ***")
        findings["status"] = "SKIP_BAD_DENOM"
        return findings
    
    zero_sij = []
    print(f"  Mandelstams mod {p}:")
    for i in range(1, N+1):
        for j in range(i+1, N+1):
            val = mandelstams[(i, j)]
            val_mod_p = F(val)
            if val_mod_p == 0:
                zero_sij.append((i, j, val))
                print(f"    s_({i},{j}) = {val} ≡ 0 (mod {p})  *** ZERO ***")
            else:
                print(f"    s_({i},{j}) = {val} ≡ {val_mod_p} (mod {p})")
    
    if not zero_sij:
        print(f"  ✓ No s_ij vanishes mod {p}")
    else:
        print(f"\n  *** {len(zero_sij)} Mandelstam invariants vanish mod {p} ***")
        findings["issues"].append({"type": "soft_collinear", "zeros": [(i, j) for (i, j, v) in zero_sij]})
    
    # Check for decoupled particles
    for a in range(1, N+1):
        all_zero = True
        for b in range(1, N+1):
            if a != b:
                key = (min(a, b), max(a, b))
                if F(mandelstams[key]) != 0:
                    all_zero = False
                    break
        if all_zero:
            print(f"  *** PARTICLE {a} IS DECOUPLED (all s_{a}* = 0 mod {p}) ***")
            findings["issues"].append({"type": "decoupled", "particle": a})
    
    # =========================================================================
    # CHECK 3: Check epsilon degeneracy
    # =========================================================================
    print("\n[CHECK 3] Epsilon degeneracy")
    print("-" * 40)
    
    eps_mod_p = F(eps)
    print(f"  epsilon(p1,p2,p3,p4) = {eps}")
    print(f"  epsilon mod {p} = {eps_mod_p}")
    
    if eps_mod_p == 0:
        print(f"  *** EPSILON VANISHES MOD {p} - BAD REDUCTION ***")
        findings["issues"].append({"type": "epsilon_zero"})
    else:
        print(f"  ✓ Epsilon non-zero mod {p}")
    
    # =========================================================================
    # Solve CHY and check solutions
    # =========================================================================
    print("\n[CHECK 4-5] Solving CHY equations")
    print("-" * 40)
    
    R, polys = build_cleared_polynomials(F, mandelstams, GAUGE)
    I_sat, solutions = saturate_and_get_solutions(R, polys, GAUGE)
    
    N_p = len(solutions)
    print(f"  Number of solutions N_p = {N_p}")
    
    if N_p == 0:
        print(f"  ✓ N_p = 0 (expected for inert prime)")
        findings["N_p"] = 0
        findings["status"] = "PASS"
    else:
        print(f"  *** N_p = {N_p} != 0 - FAILURE ***")
        findings["N_p"] = N_p
        findings["status"] = "FAIL"
        
        # Examine each solution
        for sol_idx, sol in enumerate(solutions):
            print(f"\n  Solution {sol_idx + 1}:")
            sigma_vals = {
                1: F(sigma1), 2: F(sigma2), 3: F(sigma3),
                4: sol[R.gen(0)], 5: sol[R.gen(1)], 6: sol[R.gen(2)], 7: sol[R.gen(3)]
            }
            
            print(f"    σ = ({sigma_vals[1]}, {sigma_vals[2]}, {sigma_vals[3]}, {sigma_vals[4]}, {sigma_vals[5]}, {sigma_vals[6]}, {sigma_vals[7]})")
            
            # Check for collision (should be excluded by saturation, but verify)
            collision_found = False
            for i in range(1, 8):
                for j in range(i+1, 8):
                    if sigma_vals[i] == sigma_vals[j]:
                        print(f"    *** COLLISION: σ_{i} = σ_{j} = {sigma_vals[i]} ***")
                        collision_found = True
                        findings["issues"].append({"type": "collision", "solution": sol_idx, "i": i, "j": j})
            
            if not collision_found:
                print(f"    ✓ No collisions in this solution")
            
            # Check Jacobian
            jac_det = compute_chy_jacobian(mandelstams, sigma_vals, p)
            if jac_det is None:
                print(f"    *** Jacobian computation failed (collision) ***")
            elif jac_det == 0:
                print(f"    *** JACOBIAN det = 0 - SINGULAR SOLUTION ***")
                findings["issues"].append({"type": "singular_jacobian", "solution": sol_idx})
            else:
                print(f"    Jacobian det = {jac_det} (non-zero)")
    
    return findings


def main():
    seeds = load_D1_seeds()
    
    print(f"Loaded {len(seeds)} seeds from D1")
    print(f"Investigating prime p = {P}")
    print(f"Gauge: σ₁={GAUGE[0]}, σ₂={GAUGE[1]}, σ₃={GAUGE[2]}")
    
    all_findings = []
    failures = []
    
    for seed_id, seed_value in seeds:
        findings = investigate_seed(seed_id, seed_value, P)
        all_findings.append(findings)
        
        if findings.get("status") == "FAIL":
            failures.append(findings)
    
    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    
    total_tests = len(all_findings)
    passed = sum(1 for f in all_findings if f.get("status") == "PASS")
    failed = len(failures)
    
    print(f"Total tests: {total_tests}")
    print(f"Passed (N_p = 0): {passed}")
    print(f"Failed (N_p > 0): {failed}")
    
    if failures:
        print("\nFailure details:")
        for f in failures:
            print(f"\n  Seed {f['seed_id']} (seed_value={f['seed_value']}):")
            print(f"    N_p = {f['N_p']}")
            if f.get("issues"):
                for issue in f["issues"]:
                    print(f"    Issue: {issue}")
    
    # Diagnosis
    print("\n" + "=" * 70)
    print("DIAGNOSIS")
    print("=" * 70)
    
    # Check if any common pattern
    all_issues = []
    for f in failures:
        all_issues.extend(f.get("issues", []))
    
    issue_types = {}
    for issue in all_issues:
        t = issue.get("type", "unknown")
        issue_types[t] = issue_types.get(t, 0) + 1
    
    if issue_types:
        print("Issue type frequencies:")
        for t, count in sorted(issue_types.items(), key=lambda x: -x[1]):
            print(f"  {t}: {count}")
    
    # Check if all failures have s_ij = 0
    soft_collinear_failures = [f for f in failures if any(i.get("type") == "soft_collinear" for i in f.get("issues", []))]
    if soft_collinear_failures:
        print(f"\n{len(soft_collinear_failures)} failures have s_ij ≡ 0 (mod {P})")
        print("This indicates soft/collinear degeneracy - a valid reason for bad reduction!")
    
    # Save results
    out_path = os.path.normpath(os.path.join(_BASE_DIR, "..", "logs", "p31_investigation.json"))
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    
    def json_default(obj):
        try:
            return int(obj)
        except:
            return str(obj)
    
    with open(out_path, "w") as f:
        json.dump({
            "timestamp": str(datetime.now()),
            "prime": P,
            "gauge": list(GAUGE),
            "total_tests": total_tests,
            "passed": passed,
            "failed": failed,
            "failures": failures
        }, f, indent=2, default=json_default)
    
    print(f"\nResults saved to: {out_path}")


if __name__ == "__main__":
    main()
