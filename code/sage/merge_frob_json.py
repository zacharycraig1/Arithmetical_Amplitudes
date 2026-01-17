import glob
import json


def main():
    files = sorted(glob.glob("../../data/frob_chunk_*.json"))
    if not files:
        raise SystemExit("no files matched ../../data/frob_chunk_*.json")

    out = {"records": []}
    for fn in files:
        j = json.load(open(fn))
        out["records"].extend(j.get("records", []))

    first = json.load(open(files[0]))
    for k in first:
        if k != "records":
            out[k] = first[k]

    json.dump(out, open("../../data/frob_merged_seed0.json", "w"), indent=2)
    print("merged", len(out["records"]), "records")


if __name__ == "__main__":
    main()
