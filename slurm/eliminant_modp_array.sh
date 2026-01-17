#!/usr/bin/env bash
#SBATCH -J elim
#SBATCH -N 1
#SBATCH -c 32
#SBATCH --mem=0
#SBATCH -t 12:00:00
#SBATCH -o logs/elim_%A_%a.out
#SBATCH -e logs/elim_%A_%a.err

module load sage || true

SEED_VALUE=$1
P=$2

cd code/sage

sage -python eliminant_modp_dump.sage \
  --seed_value ${SEED_VALUE} \
  --p ${P}
