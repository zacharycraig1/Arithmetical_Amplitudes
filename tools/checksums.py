import hashlib
import os
import sys


def parse_sha256sums(path):
    entries = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            digest, filename = parts[0], parts[-1]
            entries.append((digest, filename))
    return entries


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    if len(sys.argv) != 2:
        print("Usage: python tools/checksums.py data/SHA256SUMS")
        return 2

    sums_path = sys.argv[1]
    base_dir = os.path.dirname(os.path.abspath(sums_path))
    entries = parse_sha256sums(sums_path)
    if not entries:
        print("No checksum entries found.")
        return 1

    failures = 0
    for expected, filename in entries:
        target = os.path.join(base_dir, filename)
        if not os.path.exists(target):
            print(f"FAILED: {filename} (missing)")
            failures += 1
            continue
        actual = sha256_file(target)
        if actual.lower() != expected.lower():
            print(f"FAILED: {filename}")
            print(f"  expected: {expected}")
            print(f"  actual:   {actual}")
            failures += 1
        else:
            print(f"OK: {filename}")

    if failures:
        print(f"Checksum failures: {failures}")
        return 1

    print("All checksums verified.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
