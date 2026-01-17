import collections
import json


def main():
    with open("logs/n7_full.json", "r") as f:
        data = json.load(f)

    by_prime = collections.defaultdict(list)
    for t in data.get("inert_tests", []):
        by_prime[t["p"]].append(t["N_p"])

    print("\\begin{tabular}{rccc}")
    print("\\toprule")
    print("p & tests & failures & max $N_p$ \\\\")
    print("\\midrule")
    for p in sorted(by_prime):
        arr = by_prime[p]
        tests = len(arr)
        fails = sum(1 for x in arr if x != 0)
        mx = max(arr) if arr else "-"
        print(f"{p} & {tests} & {fails} & {mx} \\\\")
    print("\\bottomrule")
    print("\\end{tabular}")


if __name__ == "__main__":
    main()
