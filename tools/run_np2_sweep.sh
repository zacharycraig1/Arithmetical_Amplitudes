#!/usr/bin/env bash
set -euo pipefail

# Example sweep: inert primes x seeds
PRIMES=(7 11 19 23)
SEEDS=(0 1 2 3 4 5 6 7 8 9)

mkdir -p logs/np2

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export SAGE_NUM_THREADS=1

MAX_JOBS=${MAX_JOBS:-24}

jobcount() { jobs -pr | wc -l; }

for p in "${PRIMES[@]}"; do
  for s in "${SEEDS[@]}"; do
    while [ "$(jobcount)" -ge "$MAX_JOBS" ]; do
      sleep 0.2
    done
    echo "RUN seed=$s prime=$p"
    docker run --rm -v "$PWD:/home/sage/work" -w /home/sage/work sagemath/sagemath:10.4 \
      sage referee_checks/verify_n7_np2_extension.sage --seed "$s" --prime "$p" \
      > "logs/np2/np2_seed${s}_p${p}.txt" 2>&1 &
  done
done

wait
echo "DONE: sweep complete"
