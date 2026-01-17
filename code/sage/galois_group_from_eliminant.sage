from sage.all import *
import argparse


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--poly", type=str, default="../../data/eliminant_Z_seed42.sobj")
    args, _ = parser.parse_known_args()

    f = load(args.poly)
    R.<x> = PolynomialRing(QQ)
    f = R(f)

    print("degree", f.degree())
    try:
        print("disc factorization (partial):", factor(f.discriminant()))
    except Exception as e:
        print("disc factorization failed:", e)

    G = f.galois_group(pari_group=True)
    print("Galois group:", G)
    print("order:", G.order())

    try:
        print("structure:", G.structure_description())
    except Exception:
        pass
