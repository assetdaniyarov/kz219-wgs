#!/usr/bin/env bash
# Step 2. Joint genotyping of the 219 gVCFs into one multi-sample VCF,
# DRAGEN v4.3.6 iterative gVCF genotyper.
# Usage: bash 02_dragen_joint_genotyping.sh <gvcf_list.txt>   (one gVCF path per line)
set -euo pipefail

GVCF_LIST=$1
REF_FASTA=${REF_FASTA:?set REF_FASTA to hg19.fasta}
DRAGEN=${DRAGEN:-/opt/dragen/4.3.6/bin/dragen}
OUT=${OUT:-msvcf}
TMP=${TMP_DIR:-tmp/msvcf}
mkdir -p "$OUT" "$TMP"

"$DRAGEN" -f \
  --enable-gvcf-genotyper-iterative=true \
  --gvcfs-to-msvcf=true \
  --variant-list "$GVCF_LIST" \
  --ht-reference "$REF_FASTA" \
  --intermediate-results-dir "$TMP" \
  --output-directory "$OUT" \
  --output-file-prefix kz_219_hg19_dragen \
  --gg-msvcf-info-fields=All
