# usage:
#   sage tools/filter_transitive24_by_cycletypes.sage data/frob_cycletypes_seed0_p5_43.json

import json
import sys
from sage.all import *


def cycle_partition_of_perm(g):
    cyc = g.cycle_type()
    return Partition(sorted(cyc, reverse=True))


def main():
    fn = sys.argv[1]
    data = json.load(open(fn))
    obs_parts = sorted({tuple(r["cycle_degrees"]) for r in data["records"]})

    print("observed partitions:")
    for p in obs_parts:
        print(list(p))

    req = [Partition(list(p)) for p in obs_parts]

    candidates = []
    TG = TransitiveGroups(24)
    print("transitive groups degree 24:", len(TG))

    for k in range(1, len(TG) + 1):
        G = TransitiveGroup(24, k)
        try:
            if G.order() <= 200000:
                elems = list(G)
            else:
                elems = [G.random_element() for _ in range(5000)]
        except Exception:
            elems = [G.random_element() for _ in range(5000)]

        parts = {cycle_partition_of_perm(g) for g in elems}

        ok = True
        for r in req:
            if r not in parts:
                ok = False
                break
        if ok:
            candidates.append(k)
            print("candidate:", k, "order=", G.order())

    print("\nFINAL candidates:", candidates)


if __name__ == "__main__":
    main()
