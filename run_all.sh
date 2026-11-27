#!/usr/bin/env bash
# Runs the six notebooks in order from the repository root.
# Usage: bash run_all.sh            (all notebooks)
#        bash run_all.sh 05 06      (selected notebooks)
set -euo pipefail
cd "$(dirname "$0")"

export KZ_DRAGEN_DIR="${KZ_DRAGEN_DIR:-data/dragen}"
export KZ_PROJECT="${KZ_PROJECT:-data/popgen}"
export KZ_CALLSET="${KZ_CALLSET:-data/kz_219_hg19_dragen.annotated.header.txt}"
export KZ_SEX_FILE="${KZ_SEX_FILE:-data/kz_219_tribes.tsv}"
export KZ_WORK_DIR="${KZ_WORK_DIR:-work/03_af}"
export KZ_KNOWN_FILE="${KZ_KNOWN_FILE:-work/03_af/kz_219_hg19_AF_from_GT_known.txt}"
export KZ_CHRM_DIR="${KZ_CHRM_DIR:-data/chrM}"
export KZ_CHRY_DIR="${KZ_CHRY_DIR:-data/chrY}"

declare -A OUT=( [01]=results/01_qc [02]=results/02_popgen [03]=results/03_af
                 [04]=results/04_biomed [05]=results/05_chrM [06]=results/06_chrY )

STEPS=("$@"); [ ${#STEPS[@]} -eq 0 ] && STEPS=(01 02 03 04 05 06)
for s in "${STEPS[@]}"; do
  nb=$(ls ${s}_*.ipynb)
  echo ">>> $nb"
  KZ_OUT_DIR="${OUT[$s]}" KZ_OUT="${OUT[$s]}" \
    jupyter nbconvert --to notebook --execute --output-dir executed \
    --ExecutePreprocessor.timeout=-1 "$nb"
done
echo "done: outputs in results/, executed notebooks in executed/"
