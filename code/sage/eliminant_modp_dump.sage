"""Compute monic eliminant f_p(x) mod p and dump coefficients.

Output format per prime:
{
  "p": 101,
  "deg": 24,
  "coeffs": [c0,c1,...,c24]  # ascending order mod p
}
"""

from sage.all import *
import argparse
import json
import os

load("kinematics.sage")
load("frob_cycletypes_eliminant.sage")


def dump_for_prime(seed_value, p, outdir="../../data/eliminant_modp", gauge=(0, 1, -1)):
    os.makedirs(outdir, exist_ok=True)

    mand = mandelstams_from_seed(seed_value, 7)
    R, I = build_chy_ideal_over_fp(mand, p, gauge=gauge)
    f = eliminant_in_s4(R, I)

    if f is None or f.degree() != 24:
        return None
    if not is_squarefree_poly(f):
        return None

    f = f.monic()
    coeffs = [int(c) for c in f.coefficients(sparse=False)]

    rec = {"p": int(p), "deg": int(f.degree()), "coeffs": coeffs}
    fn = f"{outdir}/seed{seed_value}_p{p}.json"
    with open(fn, "w") as f_out:
        json.dump(rec, f_out, indent=2)
    return fn


def parse_gauge(gauge_str):
    parts = [int(x.strip()) for x in gauge_str.split(",")]
    if len(parts) != 3:
        raise ValueError("gauge must have 3 integers, e.g. 0,1,-1")
    return tuple(parts)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed_value", type=int, default=42)
    parser.add_argument("--p_min", type=int, default=None)
    parser.add_argument("--p_max", type=int, default=None)
    parser.add_argument("--outdir", type=str, default="../../data/eliminant_modp")
    parser.add_argument("--gauge", type=str, default="0,1,-1")
    args, _ = parser.parse_known_args()

    if args.p_min is None or args.p_max is None:
        raise ValueError("--p_min and --p_max are required")

    gauge = parse_gauge(args.gauge)
    for p in prime_range(args.p_min, args.p_max + 1):
        fn = dump_for_prime(args.seed_value, p, outdir=args.outdir, gauge=gauge)
        print("p", p, "->", fn)


# Only run main when executed directly
import sys
if __name__ == "__main__" and "eliminant_modp_dump" in sys.argv[0]:
    main()
