import json
from collections import Counter


def main():
    fn = "data/frob_cycletypes_seed0_p5_43.json"
    data = json.load(open(fn))
    recs = data.get("records", [])

    print("records:", len(recs))

    cnt = Counter(tuple(r["cycle_degrees"]) for r in recs)
    print("unique partitions:", len(cnt))
    for part, n in cnt.most_common(20):
        print(n, list(part))

    inert = [r for r in recs if r.get("p_mod4") == 3]
    split = [r for r in recs if r.get("p_mod4") == 1]
    print("inert count:", len(inert), "split count:", len(split))


if __name__ == "__main__":
    main()
