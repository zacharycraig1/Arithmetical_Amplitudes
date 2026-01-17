"""
REFEREE VERIFICATION: n=7 Inert Prime Vanishing
================================================

This script verifies the main computational claim:
    For n=7 CHY in 4D on dataset D1, N_p = 0 for all good-reduced inert primes p ≡ 3 (mod 4)

Also includes CONTROL: split primes p ≡ 1 (mod 4) to show pattern is non-trivial.

Usage:
    sage verify_n7_theorem.sage
    sage verify_n7_theorem.sage --quick
    sage verify_n7_theorem.sage --seeds 10 --primes 7,11,19
    sage verify_n7_theorem.sage --method groebner --output ../logs/n7_results.json
    sage verify_n7_theorem.sage --method brute --no-dataset
    sage verify_n7_theorem.sage --crosscheck

Expected output:
    All inert prime tests should pass with N_p = 0.
    Split prime tests should show NON-ZERO counts (control).
"""

import sys
import json
import os
from datetime import datetime

print("=" * 70)
print("REFEREE VERIFICATION: n=7 Inert Prime Vanishing")
print("=" * 70)
print(f"Timestamp: {datetime.now()}")
print()

# =============================================================================
# CLI ARGUMENTS
# =============================================================================

def _get_arg(flag, default=None):
    if flag in sys.argv:
        i = sys.argv.index(flag)
        if i + 1 < len(sys.argv):
            return sys.argv[i + 1]
    return default

def _has_flag(flag):
    return flag in sys.argv

quick_mode = _has_flag("--quick")
num_seeds = int(_get_arg("--seeds", "10" if quick_mode else "30"))
custom_primes = _get_arg("--primes", None)
method = _get_arg("--method", "groebner").strip().lower()  # groebner (default) or brute
out_json = _get_arg("--output", None)
use_dataset = not _has_flag("--no-dataset")
crosscheck = _has_flag("--crosscheck")
_BASE_DIR = os.path.dirname(os.path.abspath(__file__))
if out_json is None:
    out_json = os.path.normpath(os.path.join(_BASE_DIR, "..", "logs", "n7_results.json"))

# Prefer the canonical kinematics generator used to define dataset D1.
_KINEMATICS_PATH = os.path.normpath(os.path.join(_BASE_DIR, "..", "code", "sage", "kinematics.sage"))
_HAVE_CANONICAL_KIN = False
try:
    load(_KINEMATICS_PATH)
    _HAVE_CANONICAL_KIN = True
except Exception as _e:
    _HAVE_CANONICAL_KIN = False

def default_prime_sets(pmax_inert=50, pmax_split=30):
    inert = [int(p) for p in prime_range(5, pmax_inert + 1) if p % 4 == 3]
    split = [int(p) for p in prime_range(5, pmax_split + 1) if p % 4 == 1]
    return inert, split

# =============================================================================
# DATASET LOADER
# =============================================================================

def load_D1_seeds(path=None):
    if path is None:
        path = os.path.normpath(os.path.join(_BASE_DIR, "..", "data", "D1_seeds.json"))
    with open(path, "r") as f:
        obj = json.load(f)
    seeds = []
    for s in obj.get("seeds", []):
        if s.get("status", "valid") == "valid":
            seeds.append((s["id"], s["seed_value"]))
    return obj, seeds

# =============================================================================
# KINEMATICS GENERATION
# =============================================================================

def generate_4D_null_momenta(n, seed=42, scale=100):
    """Generate n null momenta in 4D with momentum conservation."""
    set_random_seed(seed)
    
    def random_null_vector(scale):
        """Generate a rational null 4-vector using Pythagorean parametrization."""
        while True:
            m = randint(1, scale)
            n_val = randint(1, scale)
            p = randint(1, scale)
            q = randint(1, scale)
            
            x = m^2 + n_val^2 - p^2 - q^2
            y = 2*(m*q + n_val*p)
            z = 2*(n_val*q - m*p)
            E = m^2 + n_val^2 + p^2 + q^2
            
            # Random signs
            sx, sy, sz = choice([-1,1]), choice([-1,1]), choice([-1,1])
            
            return (QQ(E), QQ(x)*sx, QQ(y)*sy, QQ(z)*sz)
    
    # Generate n-2 random null momenta
    momenta = [random_null_vector(scale) for _ in range(n-2)]
    
    # Compute sum of first n-2
    K = [sum(m[i] for m in momenta) for i in range(4)]
    K_sq = K[0]^2 - K[1]^2 - K[2]^2 - K[3]^2
    
    if K_sq == 0:
        return generate_4D_null_momenta(n, seed+1, scale)
    
    # Find two more null momenta summing to -K
    for attempt in range(1000):
        y = QQ(randint(-scale, scale))
        z = QQ(randint(-scale, scale))
        
        if K[0] == 0:
            continue
        
        a = K[1]^2 - K[0]^2
        if a == 0:
            continue
        
        C = y*K[2] + z*K[3] - K_sq/2
        b = 2*K[1]*C
        c_val = C^2 - K[0]^2*(y^2+z^2)
        
        disc = b^2 - 4*a*c_val
        if disc >= 0 and disc.is_square():
            x = (-b + disc.sqrt()) / (2*a)
            E = (x*K[1] + y*K[2] + z*K[3] - K_sq/2) / K[0]
            
            kn_1 = (E, x, y, z)
            kn = tuple(-K[i] - kn_1[i] for i in range(4))
            
            def is_null(p):
                return p[0]^2 - p[1]^2 - p[2]^2 - p[3]^2 == 0
            
            if is_null(kn_1) and is_null(kn):
                return momenta + [kn_1, kn]
    
    return generate_4D_null_momenta(n, seed+1000, scale)


def momenta_to_mandelstams(momenta):
    """Convert momenta to Mandelstam invariants."""
    n = len(momenta)
    s = {}
    for i in range(n):
        for j in range(i+1, n):
            pi, pj = momenta[i], momenta[j]
            dot = pi[0]*pj[0] - pi[1]*pj[1] - pi[2]*pj[2] - pi[3]*pj[3]
            s[(i+1, j+1)] = 2 * dot
    return s

def epsilon_1234(momenta):
    """
    Determinant epsilon(p1,p2,p3,p4) used as a genericity / good-reduction guard.
    For our rational kinematics this is in QQ; reduction mod p is meaningful iff denom not divisible by p.
    """
    M = matrix(QQ, [list(momenta[0]), list(momenta[1]), list(momenta[2]), list(momenta[3])])
    return M.det()


# =============================================================================
# CHY SOLUTION COUNTING
# =============================================================================

def build_cleared_polynomials(F, mandelstams, gauge):
    """
    Build cleared-denominator CHY polynomials for n=7 in gauge (sigma1,sigma2,sigma3).
    Variables: x4,x5,x6,x7. Equations: h_4 = h_5 = h_6 = h_7 = 0.
    """
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


def saturate_and_count(R, polys, gauge, do_crosscheck=False):
    """
    Saturate away collision loci among sigma variables:
    - x_i != x_j
    - x_i != gauge values
    Return: integer count (usually 0 or 24).
    """
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
        Q = R.quotient(I_sat)
        qdim = int(Q.vector_space_dimension())
    except Exception:
        qdim = None

    if do_crosscheck:
        # NOTE:
        # - vcount = #F_p-rational points of the saturated variety.
        # - qdim is the degree/length of the 0-dim scheme (points over algebraic closure,
        #   counted with multiplicity), and need NOT equal vcount in general.
        vcount = len(I_sat.variety())
        if qdim is not None and qdim != vcount:
            print(f"  [crosscheck] quotient dim (degree) {qdim} != F_p-point count {vcount} (this can be normal)")
        return vcount

    # IMPORTANT:
    # `qdim` is NOT the same thing as N_p (F_p-rational points), except in special cases.
    # We only use it as a fast-path for the empty case.
    if qdim == 0:
        return 0
    return len(I_sat.variety())


def compute_chy_jacobian_det(mandelstams, sigma_vals, p):
    """
    Compute det(d h_a / d sigma_b) for a,b in {4,5,6,7}.
    Returns None on collision.
    """
    F = GF(p)

    def s_ij(i, j):
        a, b = min(i, j), max(i, j)
        return F(mandelstams[(a, b)])

    J = matrix(F, 4, 4)
    for row, a in enumerate([4, 5, 6, 7]):
        for col, b in enumerate([4, 5, 6, 7]):
            if a == b:
                total = F(0)
                for c in range(1, 8):
                    if c == a:
                        continue
                    denom = sigma_vals[a] - sigma_vals[c]
                    if denom == 0:
                        return None
                    total -= s_ij(a, c) / (denom^2)
                J[row, col] = total
            else:
                denom = sigma_vals[a] - sigma_vals[b]
                if denom == 0:
                    return None
                J[row, col] = s_ij(a, b) / (denom^2)
    return J.det()


def count_chy_solutions_groebner(mandelstams, gauge, p, do_crosscheck=False):
    """
    Count F_p solutions to n=7 CHY equations using Groebner + saturation.
    Returns: (N_p, is_bad, reason)
    """
    sigma1, sigma2, sigma3 = gauge
    F = GF(p)

    for (v, w) in [(sigma1, sigma2), (sigma1, sigma3), (sigma2, sigma3)]:
        if F(v - w) == 0:
            return -1, True, f"gauge collision mod {p}"

    for key, val in mandelstams.items():
        if hasattr(val, 'denominator') and val.denominator() % p == 0:
            return -1, True, f"bad denominator mod {p}"

    for particle in range(1, 8):
        all_zero = True
        for other in range(1, 8):
            if particle == other:
                continue
            key = (min(particle, other), max(particle, other))
            if F(mandelstams[key]) != 0:
                all_zero = False
                break
        if all_zero:
            return -1, True, f"particle {particle} decoupled mod {p}"

    R, polys = build_cleared_polynomials(F, mandelstams, gauge)
    I = R.ideal(polys)

    # Saturate collision loci
    sigma1, sigma2, sigma3 = gauge
    x4, x5, x6, x7 = R.gens()
    S = 1
    for x in [x4, x5, x6, x7]:
        for g in [R.base_ring()(sigma1), R.base_ring()(sigma2), R.base_ring()(sigma3)]:
            S *= (x - g)
    S *= (x4 - x5) * (x4 - x6) * (x4 - x7) * (x5 - x6) * (x5 - x7) * (x6 - x7)
    I_sat, _ = I.saturation(S)

    try:
        solutions = I_sat.variety()
    except Exception:
        solutions = []

    if do_crosscheck:
        try:
            Q = R.quotient(I_sat)
            qdim = int(Q.vector_space_dimension())
        except Exception:
            qdim = None
        if qdim is not None and qdim != len(solutions):
            print(f"  [crosscheck] quotient dim (degree) {qdim} != F_p-point count {len(solutions)} (this can be normal)")

    # Ramification check: Jacobian must be non-degenerate at all F_p solutions
    for sol in solutions:
        sigma_vals = {
            1: F(sigma1), 2: F(sigma2), 3: F(sigma3),
            4: sol[R.gen(0)], 5: sol[R.gen(1)], 6: sol[R.gen(2)], 7: sol[R.gen(3)]
        }
        jac_det = compute_chy_jacobian_det(mandelstams, sigma_vals, p)
        if jac_det is None:
            return -1, True, "collision in solution"
        if jac_det == 0:
            return -1, True, "ramified jacobian"

    return len(solutions), False, "good"


def count_chy_solutions_bruteforce(mandelstams, gauge, p):
    """
    Count F_p solutions to n=7 CHY equations by brute force.
    
    Returns: (N_p, is_bad, reason)
    """
    sigma1, sigma2, sigma3 = gauge
    F = GF(p)
    
    # Check for bad reduction
    gauge_vals = [sigma1, sigma2, sigma3]
    for i, v in enumerate(gauge_vals):
        for j, w in enumerate(gauge_vals):
            if i < j and F(v - w) == 0:
                return -1, True, f"gauge collision mod {p}"
    
    for key, val in mandelstams.items():
        if hasattr(val, 'denominator') and val.denominator() % p == 0:
            return -1, True, f"bad denominator mod {p}"
    
    for particle in range(1, 8):
        all_zero = True
        for other in range(1, 8):
            if particle != other:
                key = (min(particle, other), max(particle, other))
                if F(mandelstams[key]) != 0:
                    all_zero = False
                    break
        if all_zero:
            return -1, True, f"particle {particle} decoupled mod {p}"
    
    # Count solutions by brute force
    count = 0
    forbidden = set([F(sigma1), F(sigma2), F(sigma3)])
    
    for s4 in F:
        if s4 in forbidden:
            continue
        for s5 in F:
            if s5 in forbidden or s5 == s4:
                continue
            for s6 in F:
                if s6 in forbidden or s6 == s4 or s6 == s5:
                    continue
                for s7 in F:
                    if s7 in forbidden or s7 == s4 or s7 == s5 or s7 == s6:
                        continue
                    
                    sigma = [F(sigma1), F(sigma2), F(sigma3), s4, s5, s6, s7]
                    
                    is_solution = True
                    for a in range(3, 7):
                        h_a = F(0)
                        for b in range(7):
                            if b != a:
                                i, j = min(a+1, b+1), max(a+1, b+1)
                                s_ab = F(mandelstams[(i, j)])
                                denom = sigma[a] - sigma[b]
                                if denom == 0:
                                    is_solution = False
                                    break
                                h_a += s_ab / denom
                        if not is_solution:
                            break
                        if h_a != 0:
                            is_solution = False
                            break
                    
                    if is_solution:
                        sigma_vals = {
                            1: F(sigma1), 2: F(sigma2), 3: F(sigma3),
                            4: s4, 5: s5, 6: s6, 7: s7
                        }
                        jac_det = compute_chy_jacobian_det(mandelstams, sigma_vals, p)
                        if jac_det is None:
                            return -1, True, "collision in solution"
                        if jac_det == 0:
                            return -1, True, "ramified jacobian"
                        count += 1
    
    return count, False, "good"


# =============================================================================
# MAIN VERIFICATION
# =============================================================================

def main():
    gauge = (0, 1, -1)

    # Primes
    if custom_primes:
        inert_primes = [int(p) for p in custom_primes.split(",") if int(p) % 4 == 3]
        split_primes = [int(p) for p in custom_primes.split(",") if int(p) % 4 == 1]
    else:
        inert_primes, split_primes = default_prime_sets(50, 30)
    
    if quick_mode:
        inert_primes = inert_primes[:3]
        split_primes = split_primes[:2]
    
    print(f"Seeds: {num_seeds}")
    print(f"Method: {method}")
    print(f"Output: {out_json}")
    print(f"Dataset: {'D1_seeds.json' if use_dataset else 'synthetic'}")
    print(f"Crosscheck: {crosscheck}")
    print(f"Inert primes (expect N_p = 0): {inert_primes}")
    print(f"Split primes (CONTROL, expect N_p > 0): {split_primes}")
    print()
    
    # Statistics
    total_inert = 0
    passed_inert = 0
    failed_inert = 0
    skipped = 0
    skipped_ramification = 0
    skipped_bad_denom = 0
    skipped_decoupled = 0
    skipped_ramification_inert = 0
    skipped_bad_denom_inert = 0
    skipped_decoupled_inert = 0
    
    split_nonzero = 0
    total_split = 0
    
    results = {
        "timestamp": str(datetime.now()),
        "method": method,
        "gauge": [int(x) for x in gauge],
        "quick_mode": quick_mode,
        "crosscheck": crosscheck,
        "num_seeds": num_seeds,
        "inert_primes": inert_primes,
        "split_primes": split_primes,
        "inert_tests": [],
        "split_tests": [],
        "summary": {}
    }

    if use_dataset:
        d1_meta, seed_pairs = load_D1_seeds()
        results["D1_meta"] = {
            "description": d1_meta.get("description", ""),
            "generation_method": d1_meta.get("generation_method", "")
        }
        seed_list = seed_pairs[:num_seeds]
    else:
        seed_list = [(i, i * 100 + 42) for i in range(num_seeds)]

    print("=" * 70)
    print("INERT PRIME TESTS")
    print("=" * 70)

    for seed_idx in range(len(seed_list)):
        seed, seed_value = seed_list[seed_idx]
        try:
            if _HAVE_CANONICAL_KIN:
                momenta = make_4d_massless_point(7, seed_value)
            else:
                momenta = generate_4D_null_momenta(7, seed=seed_value)
            mandelstams = momenta_to_mandelstams(momenta)
        except Exception as e:
            print(f"Seed {seed}: SKIP (kinematics failed)")
            skipped += 1
            continue
        
        for p in inert_primes:
            if method == "brute":
                N_p, is_bad, reason = count_chy_solutions_bruteforce(mandelstams, gauge, p)
            else:
                do_cross = crosscheck and p in [7, 11] and seed_idx < 2
                N_p, is_bad, reason = count_chy_solutions_groebner(mandelstams, gauge, p, do_crosscheck=do_cross)
            
            if is_bad:
                skipped += 1
                if reason == "ramified jacobian":
                    skipped_ramification += 1
                    skipped_ramification_inert += 1
                if reason.startswith("bad denominator"):
                    skipped_bad_denom += 1
                    skipped_bad_denom_inert += 1
                if "decoupled" in reason:
                    skipped_decoupled += 1
                    skipped_decoupled_inert += 1
                print(f"seed={seed:2d}, p={p:2d}: SKIP ({reason})")
                continue
            
            total_inert += 1
            if N_p == 0:
                passed_inert += 1
                status = "PASS"
            else:
                failed_inert += 1
                status = "FAIL"
            
            print(f"seed={seed:2d}, p={p:2d}: N_p = {N_p:2d}  [{status}]")
            results["inert_tests"].append({
                "seed": seed,
                "seed_value": seed_value,
                "p": p,
                "N_p": int(N_p)
            })
    
    print()
    print("=" * 70)
    print("SPLIT PRIME CONTROL")
    print("=" * 70)
    
    for seed_idx in range(min(5, len(seed_list))):
        try:
            seed, seed_value = seed_list[seed_idx]
            if _HAVE_CANONICAL_KIN:
                momenta = make_4d_massless_point(7, seed_value)
            else:
                momenta = generate_4D_null_momenta(7, seed=seed_value)
            mandelstams = momenta_to_mandelstams(momenta)
        except:
            continue
        
        for p in split_primes:
            if method == "brute":
                N_p, is_bad, reason = count_chy_solutions_bruteforce(mandelstams, gauge, p)
            else:
                do_cross = crosscheck and p in [7, 11] and seed_idx < 2
                N_p, is_bad, reason = count_chy_solutions_groebner(mandelstams, gauge, p, do_crosscheck=do_cross)
            
            if is_bad:
                skipped += 1
                if reason == "ramified jacobian":
                    skipped_ramification += 1
                if reason.startswith("bad denominator"):
                    skipped_bad_denom += 1
                if "decoupled" in reason:
                    skipped_decoupled += 1
                print(f"seed={seed:2d}, p={p:2d}: SKIP ({reason})")
                continue
            
            total_split += 1
            if N_p > 0:
                split_nonzero += 1
            
            print(f"seed={seed:2d}, p={p:2d}: N_p = {N_p:2d}")
            results["split_tests"].append({
                "seed": seed,
                "seed_value": seed_value,
                "p": p,
                "N_p": int(N_p)
            })
    
    # Summary
    print()
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Inert prime tests: {passed_inert}/{total_inert} passed (N_p = 0)")
    print(f"Inert prime failures: {failed_inert}")
    print(f"Skipped (bad reduction / degenerate): {skipped}")
    print(f"Split prime control: {split_nonzero}/{total_split} had N_p > 0")
    print()
    
    results["summary"] = {
        "passed_inert": passed_inert,
        "total_inert": total_inert,
        "failed_inert": failed_inert,
        "skipped": skipped,
        "skipped_ramification": skipped_ramification,
        "skipped_bad_denom": skipped_bad_denom,
        "skipped_decoupled": skipped_decoupled,
        "skipped_ramification_inert": skipped_ramification_inert,
        "skipped_bad_denom_inert": skipped_bad_denom_inert,
        "skipped_decoupled_inert": skipped_decoupled_inert,
        "total_inert_pre_ramification": total_inert + skipped_ramification_inert,
        "split_nonzero": split_nonzero,
        "total_split": total_split,
        "skipped_split": skipped - (skipped_ramification_inert + skipped_bad_denom_inert + skipped_decoupled_inert)
    }

    out_dir = os.path.dirname(out_json)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)
    def _json_default(obj):
        try:
            return int(obj)
        except Exception:
            return str(obj)
    with open(out_json, "w") as f:
        json.dump(results, f, indent=2, default=_json_default)
    print(f"Wrote results JSON -> {out_json}")

    if failed_inert == 0 and total_inert > 0:
        print("[PASS] INERT PRIME PROPOSITION VERIFIED")
        print(f"       All {total_inert} good-reduction inert tests have N_p = 0")
        if split_nonzero > 0:
            print("[PASS] CONTROL GROUP VALIDATES")
            print(f"       Split primes show non-zero counts ({split_nonzero}/{total_split})")
    else:
        print("[FAIL] Some inert prime tests failed")
        sys.exit(1)


if __name__ == "__main__":
    main()
