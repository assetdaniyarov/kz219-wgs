#!/usr/bin/env bash
# Step 5. Y-chromosome haplogroups of the 74 XY genomes, Y-LineageTracker v1.3.0
# (ISOGG 2019 tree shipped with that release, GRCh37 coordinates, `classify -b 37`).
# The chrY calls come from the same DRAGEN v3.10.16 / GRCh37 callset as the autosomes;
# the XY samples are taken from the DRAGEN ploidy call and passed as an explicit list.
# Outputs kz_219_ylt.hapresult.hg and kz_219_ylt.lineageresult.txt are the inputs of
# notebook 06 (KZ_CHRY_DIR).
# Usage: bash 05_chrY_ylt.sh <multisample.vcf.gz> <xy_samples.txt>
set -euo pipefail

VCF=$1; XY=$2
OUT=${OUT:-ylt}
BCFTOOLS=${BCFTOOLS:-bcftools}
LINEAGETRACKER=${LINEAGETRACKER:-LineageTracker}
mkdir -p "$OUT"

"$BCFTOOLS" view -S "$XY" -r chrY -Oz -o "$OUT/kz_219_chrY_xy.vcf.gz" "$VCF"
"$BCFTOOLS" index -t "$OUT/kz_219_chrY_xy.vcf.gz"

"$LINEAGETRACKER" classify \
  --vcf "$OUT/kz_219_chrY_xy.vcf.gz" \
  -b 37 \
  -o "$OUT/kz_219_ylt"
