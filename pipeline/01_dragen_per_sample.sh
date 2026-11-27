#!/usr/bin/env bash
# Step 1. Per-sample alignment and gVCF calling, DRAGEN v3.10.16, GRCh37/hg19.
# Usage: bash 01_dragen_per_sample.sh <sample_id> <R1.fq.gz> <R2.fq.gz>
set -euo pipefail

SAMPLE=$1; R1=$2; R2=$3
REF_DIR=${REF_DIR:?set REF_DIR to the DRAGEN hg19 hash table}
OUT_ROOT=${OUT_ROOT:-gvcf}
LIBRARY=${LIBRARY:-bgi_kz}

OUT=$OUT_ROOT/$SAMPLE
mkdir -p "$OUT"

dragen \
  -r "$REF_DIR" \
  -1 "$R1" -2 "$R2" \
  --output-directory "$OUT" \
  --output-file-prefix "$SAMPLE" \
  --RGID "$LIBRARY" --RGLB "$LIBRARY" --RGPL Illumina --RGSM "$SAMPLE" --RGCN NLA \
  --enable-map-align true \
  --enable-map-align-output true \
  --enable-sort true \
  --enable-duplicate-marking true \
  --enable-bam-indexing true \
  --generate-md-tags true \
  --preserve-bqsr-tags true \
  --enable-variant-caller true \
  --vc-emit-ref-confidence GVCF \
  --vc-enable-phasing true \
  --vc-enable-roh true \
  --vc-enable-af-filter true \
  --enable-cnv true \
  --enable-sv true \
  --repeat-genotype-enable true
