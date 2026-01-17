#!/usr/bin/env bash
set -euo pipefail

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export NUMEXPR_MAX_THREADS=1
export SAGE_NUM_THREADS=1
export TMPDIR=/dev/shm

make USE_DOCKER=0 checksums
make USE_DOCKER=0 verify-full
make USE_DOCKER=0 verify-big

echo "ALL DONE. Logs are in logs/"
