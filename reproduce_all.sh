#!/usr/bin/env bash
set -euo pipefail

# One-command reproduction (offline-friendly)

make docker-pull
make checksums
make verify-full

echo "DONE: paper claims reproduced. Logs in logs/"
