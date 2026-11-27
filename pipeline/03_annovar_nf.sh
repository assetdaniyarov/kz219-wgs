#!/usr/bin/env bash
# Step 3. Annotation of the multi-sample VCF with ANNOVAR v2025Mar02 via annovar-nf
# (https://github.com/assetdaniyarov/annovar-nf). The output table is the input of
# notebook 03 (KZ_CALLSET).
# Usage: bash 03_annovar_nf.sh <msvcf_dir>
set -euo pipefail

INPUT=$1
ANNOVAR_NF=${ANNOVAR_NF:-annovar-nf}          # clone of the annovar-nf repository
ANNOVAR_DB=${ANNOVAR_DB:?set ANNOVAR_DB to annovar humandb}
ANNOVAR_SRC=${ANNOVAR_SRC:?set ANNOVAR_SRC to the annovar_20250302 install}
THREADS=${THREADS:-32}

nextflow run "$ANNOVAR_NF/main.nf" \
  --input       "$INPUT" \
  --genome      hg19 \
  --annovar_db  "$ANNOVAR_DB" \
  --annovar_src "$ANNOVAR_SRC" \
  --threads     "$THREADS" \
  --outdir      "${OUTDIR:-$INPUT}"
