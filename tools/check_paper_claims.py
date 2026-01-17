import json


def main() -> None:
    with open("logs/n7_full.json", "r") as handle:
        data = json.load(handle)
    summary = data["summary"]

    assert summary["passed_inert"] == 183
    assert summary["failed_inert"] == 0
    assert summary["skipped_bad_denom_inert"] == 11
    assert summary["skipped_ramification_inert"] == 3
    assert summary["skipped_decoupled_inert"] == 13
    assert summary["total_inert_pre_ramification"] == 186

    print("OK: paper-facing inert-prime counts match logs/n7_full.json")


if __name__ == "__main__":
    main()
