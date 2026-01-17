"""Frobenius cycle types via eliminant factor degrees (n=7 CHY).

This avoids solution-labeling issues by extracting cycle types from the
factor degrees of a univariate eliminant modulo p.

Good reduction criterion used here:
- eliminant exists
- degree = 24
- squarefree mod p
"""

from sage.all import *
import argparse
import json
import os

# ------------------------------------------------------------
# 1) Kinematics (deterministic generator)
# ------------------------------------------------------------
load("kinematics.sage")


def mandelstams_from_seed(seed_value, n=7):
    momenta = make_4d_massless_point(n, seed_value)
    return momenta_to_mandelstams(momenta)


# ------------------------------------------------------------
# 2) Build CHY polynomial system over F_p
# ------------------------------------------------------------
def build_chy_ideal_over_fp(mandelstams, p, gauge=(0, 1, -1)):
    """Return (R, I) where I is the CHY ideal in R = GF(p)[s4,s5,s6,s7].

    We fix sigma1, sigma2, sigma3 to gauge values and clear denominators.
    Equations: h_a = sum_{b!=a} s_ab/(sigma_a - sigma_b) = 0 for a=4..7.
    """
    F = GF(p)

    # Lex order so that elimination gives a univariate in s4
    R = PolynomialRing(F, names=["s7", "s6", "s5", "s4"], order="lex")
    s7, s6, s5, s4 = R.gens()

    # Gauge-fixed sigmas
    s1, s2, s3 = [F(g) for g in gauge]
    sig = {1: s1, 2: s2, 3: s3, 4: s4, 5: s5, 6: s6, 7: s7}

    # Precompute pairwise Mandelstams mod p
    def s_ab(a, b):
        i, j = min(a, b), max(a, b)
        return F(mandelstams[(i, j)])

    polys = []
    for a in [4, 5, 6, 7]:
        # denom_a = ∏_{b!=a} (σ_a - σ_b)
        denom_a = R(1)
        for b in [1, 2, 3, 4, 5, 6, 7]:
            if b == a:
                continue
            denom_a *= (sig[a] - sig[b])

        # numerator after clearing denominators
        num = R(0)
        for b in [1, 2, 3, 4, 5, 6, 7]:
            if b == a:
                continue
            term = s_ab(a, b)
            num += term * (denom_a // (sig[a] - sig[b]))

        polys.append(num)

    I = Ideal(polys)
    return R, I


# ------------------------------------------------------------
# 3) Extract eliminant f_p(s4)
# ------------------------------------------------------------
def eliminant_in_s4(R, I):
    """Return univariate eliminant in s4 by eliminating s7,s6,s5.

    We use elimination ideal: I ∩ GF(p)[s4].
    In lex order (s7 > s6 > s5 > s4), elimination is stable.
    """
    s7, s6, s5, s4 = R.gens()

    # Eliminate s7,s6,s5; keep s4
    J = I.elimination_ideal([s7, s6, s5])
    gens = J.gens()

    # Search for a nonzero univariate polynomial in s4
    univs = []
    for g in gens:
        if set(g.variables()) <= set([s4]) and g != 0:
            univs.append(g)

    if len(univs) == 0:
        return None

    # Use the lowest-degree nonzero univariate
    f = min(univs, key=lambda h: h.degree())

    # Normalize to monic if possible
    try:
        f = f.monic()
    except Exception:
        pass

    return f


# ------------------------------------------------------------
# 4) Good reduction + cycle type extraction
# ------------------------------------------------------------
def cycle_type_from_factor_degrees(f):
    """Return sorted list of irreducible factor degrees (descending)."""
    facs = f.factor()
    degs = []
    for g, e in facs:
        degs += [g.degree()] * e
    degs.sort(reverse=True)
    return degs


def is_squarefree_poly(f):
    """Squarefree test via gcd(f, f')."""
    if f.degree() <= 0:
        return False
    fp = f.derivative()
    return gcd(f, fp) == 1


def classify_cycle_degs_D24_rotation(degs):
    """Return True if degrees match one of the D24 rotation templates."""
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


# ------------------------------------------------------------
# 5) Main driver
# ------------------------------------------------------------
def run(
    seed_id=0,
    seed_value=42,
    p_min=5,
    p_max=500,
    outfile="../../data/frob_cycletypes_seed0.json",
    gauge=(0, 1, -1),
):
    mand = mandelstams_from_seed(seed_value, 7)

    results = {
        "seed_id": seed_id,
        "seed_value": seed_value,
        "gauge": list(gauge),
        "p_range": [p_min, p_max],
        "records": [],
    }

    bad_log = []

    for p in prime_range(p_min, p_max + 1):
        if p < 5:
            continue

        try:
            R, I = build_chy_ideal_over_fp(mand, p, gauge=gauge)
            f = eliminant_in_s4(R, I)
            if f is None:
                bad_log.append((p, "no eliminant"))
                continue
            if f.degree() != 24:
                bad_log.append((p, f"wrong degree {f.degree()}"))
                continue
            if not is_squarefree_poly(f):
                bad_log.append((p, "not squarefree"))
                continue
            # Stronger bad-reduction filters (catch hidden multiplicities)
            try:
                f_sf = f.squarefree_part()
                if f_sf.degree() != 24:
                    bad_log.append((p, f"squarefree_part degree {f_sf.degree()}"))
                    continue
            except Exception:
                pass
            try:
                if f.discriminant() == 0:
                    bad_log.append((p, "disc=0"))
                    continue
            except Exception:
                pass

            degs = cycle_type_from_factor_degrees(f)
            # Convert Sage Integers to Python int for JSON serialization
            degs = [int(d) for d in degs]

            rec = {
                "p": int(p),
                "p_mod4": int(p % 4),
                "cycle_degrees": degs,
                "D24_rotation_template": classify_cycle_degs_D24_rotation(degs),
            }
            results["records"].append(rec)
        except Exception as e:
            bad_log.append((p, "exception: " + str(e)))
            continue

    os.makedirs(os.path.dirname(outfile), exist_ok=True)
    with open(outfile, "w") as f_out:
        json.dump(results, f_out, indent=2)

    badfile = outfile.replace(".json", "_bad_reduction.txt")
    with open(badfile, "w") as f_bad:
        for p, reason in bad_log:
            f_bad.write(f"p={p} :: {reason}\n")

    print("Wrote", outfile)
    print("Wrote", badfile)


def parse_gauge(gauge_str):
    parts = [int(x.strip()) for x in gauge_str.split(",")]
    if len(parts) != 3:
        raise ValueError("gauge must have 3 integers, e.g. 0,1,-1")
    return tuple(parts)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed_id", type=int, default=0)
    parser.add_argument("--seed_value", type=int, default=42)
    parser.add_argument("--p_min", type=int, default=5)
    parser.add_argument("--p_max", type=int, default=500)
    parser.add_argument("--outfile", type=str, default="../../data/frob_cycletypes_seed0.json")
    parser.add_argument("--gauge", type=str, default="0,1,-1")
    args, _ = parser.parse_known_args()

    run(
        seed_id=args.seed_id,
        seed_value=args.seed_value,
        p_min=args.p_min,
        p_max=args.p_max,
        outfile=args.outfile,
        gauge=parse_gauge(args.gauge),
    )


# Only run main when executed directly, not when load()'ed by Sage
import sys
if __name__ == "__main__" and "frob_cycletypes_eliminant" in sys.argv[0]:
    main()
