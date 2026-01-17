# usage:
#   sage tools/filter_by_derangement_coset.sage <k>
# where k is 24T{k}

import sys
from sage.all import *


def is_derangement(g, omega):
    return all(g(i) != i for i in omega)


def main():
    k = int(sys.argv[1])
    G = TransitiveGroup(24, k)
    print("G=24T%d order=%d" % (k, G.order()))

    omega = list(range(1, 25))
    subs = [H for H in G.subgroups() if H.order() * 2 == G.order()]
    print("index-2 subgroups:", len(subs))

    for idx, H in enumerate(subs):
        coset_rep = None
        for g in G.gens():
            if g not in H:
                coset_rep = g
                break
        if coset_rep is None:
            for _ in range(1000):
                g = G.random_element()
                if g not in H:
                    coset_rep = g
                    break
        if coset_rep is None:
            continue

        good = 0
        bad = 0
        for _ in range(2000):
            h = H.random_element()
            g = coset_rep * h
            if is_derangement(g, omega):
                good += 1
            else:
                bad += 1

        if bad <= 5:
            print("MATCH: subgroup", idx, "bad", bad, "good", good)


if __name__ == "__main__":
    main()
