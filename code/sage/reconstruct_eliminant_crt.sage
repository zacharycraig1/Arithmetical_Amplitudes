"""CRT reconstruct eliminant coefficients into Z[x].

Assumes each modp json contains monic degree-24 polynomial.
"""

from sage.all import *
import argparse
import glob
import json


def load_modp_records(pattern):
    files = sorted(glob.glob(pattern))
    recs = []
    for fn in files:
        recs.append(json.load(open(fn)))
    return recs


def crt_reconstruct_coeff_vector(recs):
    deg = recs[0]["deg"]
    ncoeff = deg + 1

    M = 1
    coeffs_modM = [0] * ncoeff

    for rec in recs:
        p = Integer(rec["p"])
        a = [Integer(x) for x in rec["coeffs"]]
        if len(a) != ncoeff:
            raise ValueError("wrong coeff length")

        new_coeffs = []
        for i in range(ncoeff):
            new_coeffs.append(crt(coeffs_modM[i], a[i], M, p))

        coeffs_modM = new_coeffs
        M *= p

    return M, coeffs_modM


def center_lift(vec, M):
    out = []
    for x in vec:
        y = Integer(x)
        if y > M // 2:
            y -= M
        out.append(y)
    return out


def build_poly_from_coeffs(coeffs):
    R.<x> = PolynomialRing(ZZ)
    f = sum(ZZ(coeffs[i]) * x**i for i in range(len(coeffs)))
    return f


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--pattern", type=str, default="../../data/eliminant_modp/seed42_p*.json")
    parser.add_argument("--out_json", type=str, default="../../data/eliminant_Z_seed42.json")
    parser.add_argument("--out_poly", type=str, default="../../data/eliminant_Z_seed42.sobj")
    args, _ = parser.parse_known_args()

    recs = load_modp_records(args.pattern)
    if len(recs) == 0:
        raise ValueError("no input files matched pattern")

    M, coeffs_modM = crt_reconstruct_coeff_vector(recs)
    coeffs = center_lift(coeffs_modM, M)
    f = build_poly_from_coeffs(coeffs)

    out = {
        "modulus_M": str(M),
        "coeffs_centered": [int(c) for c in coeffs],
    }
    json.dump(out, open(args.out_json, "w"), indent=2)
    f.save(args.out_poly)

    print("loaded", len(recs), "primes")
    print("CRT modulus bits:", M.nbits())
    print("degree:", f.degree())
    print("leading coeff:", f.leading_coefficient())
    print("saved", args.out_json, "and", args.out_poly)
