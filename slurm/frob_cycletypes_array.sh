#!/usr/bin/env bash
#SBATCH -J frob
#SBATCH -N 1
#SBATCH -c 32
#SBATCH --mem=0
#SBATCH -t 08:00:00
#SBATCH -o logs/frob_%A_%a.out
#SBATCH -e logs/frob_%A_%a.err

module load sage || true

SEED_VALUE=$1
CHUNK_MIN=$2
CHUNK_MAX=$3
OUTFILE=$4

cd code/sage

sage -python frob_cycletypes_eliminant.sage \
  --seed_value ${SEED_VALUE} \
  --p_min ${CHUNK_MIN} \
  --p_max ${CHUNK_MAX} \
  --outfile ${OUTFILE}
