#!/bin/bash
set -euo pipefail
source "$(dirname "$0")/config.sh"

# ============================================================================
# gq_profile.sh — decide KZ_GQ_MIN from the data, before launching 01.
#
# WHY THIS EXISTS. check_inputs.sh reported "GQ < 20: 75.7% of the first 165,090
# genotype calls". That number is not usable: it came from `head` on a stream
# that starts at the beginning of chr1, i.e. the first ~750 variants of the
# assembly, which sit in the telomeric/subtelomeric pile-up where essentially
# every callset has poor confidence. Taken at face value it says GQ masking is
# impossible; measured properly it may say the opposite. Neither guess is worth
# a day of compute.
#
# WHAT IT DOES. Samples PASS sites from several 200 kb windows spread across
# autosomes, away from centromeres and telomeres, and then REPLAYS THE EXACT
# FILTER CHAIN 01 uses:
#
#     mask genotypes with GQ < T  ->  drop variants with missingness > GENO_PRE
#                                 ->  measure per-sample missingness on what is
#                                     left  ->  would --mind MIND_WGS fire?
#
# It does this for several candidate thresholds at once, including T=0, which is
# the no-genotype-filter baseline. The output is a table you can read a decision
# off: the largest T whose variant loss and per-sample missingness you are
# willing to pay for.
#
# This is a sample of ~1.6 Mb, not the genome. It is representative enough to
# separate "GQ 20 is fine" from "GQ 20 deletes the panel", which is the only
# question being asked here.
#
#   bash gq_profile.sh
#   GQ_THRESHOLDS="0 10 15 20 30" bash gq_profile.sh
# ============================================================================

command -v bcftools >/dev/null || die "bcftools not found; needed for the GQ profile"
require_file "${DS_KZ_VCF}"

GQ_THRESHOLDS="${GQ_THRESHOLDS:-0 10 15 20 30}"
OUT="${RES_DIR}/provenance"; mkdir -p "${OUT}"

# An index is required: without one bcftools cannot seek and we are back to
# reading the head of chr1, which is the whole problem.
if [[ ! -s "${DS_KZ_VCF}.tbi" && ! -s "${DS_KZ_VCF}.csi" ]]; then
    die "no .tbi/.csi beside ${DS_KZ_VCF}. Run: bcftools index -t '${DS_KZ_VCF}'"
fi

# Contig naming: '1' or 'chr1'. Read it off the index rather than assuming.
PFX=""
bcftools index -s "${DS_KZ_VCF}" 2>/dev/null | head -1 | grep -q '^chr' && PFX="chr"

# 200 kb windows on eight autosomes, all in gene-containing euchromatin, none
# adjacent to a centromere or telomere (GRCh37 coordinates).
REGIONS="${GQ_REGIONS:-}"
if [[ -z "${REGIONS}" ]]; then
    for R in 1:20000000-20200000 2:50000000-50200000 4:100000000-100200000 \
             7:80000000-80200000 10:60000000-60200000 14:60000000-60200000 \
             18:40000000-40200000 20:35000000-35200000; do
        REGIONS="${REGIONS:+${REGIONS},}${PFX}${R}"
    done
fi
log "=== GQ profile ==="
log "  VCF:        ${DS_KZ_VCF}"
log "  regions:    ${REGIONS}"
log "  thresholds: ${GQ_THRESHOLDS}"
log "  replaying:  mask GQ<T -> --geno ${GENO_PRE} -> --mind ${MIND_KZ}"

# PASS only, to match KZ_VCF_FILTERS (--var-filter). One line per variant,
# one field per sample, "GT:GQ".
bcftools view -f PASS -r "${REGIONS}" -Ou "${DS_KZ_VCF}" 2>/dev/null \
  | bcftools query -f '%CHROM:%POS[\t%GT:%GQ]\n' \
  > "${WORK_DIR}/gq_sample.txt" \
  || die "bcftools failed; check the region names against 'bcftools index -s'"

NV=$(nlines "${WORK_DIR}/gq_sample.txt")
[[ "${NV}" -gt 100 ]] || die "only ${NV} PASS variants in the sampled regions"
log "  sampled: ${NV} PASS variants"

awk -v THR="${GQ_THRESHOLDS}" -v GENO_PRE="${GENO_PRE}" -v MIND="${MIND_KZ}" '
    BEGIN { FS="\t"; nt = split(THR, T, " ") }
    {
        nv++
        ns = NF - 1
        for (t = 1; t <= nt; t++) {
            miss = 0
            for (s = 2; s <= NF; s++) {
                split($s, g, ":")
                gt = g[1]; gq = g[2]
                # already missing, or masked by this threshold
                if (gt == "./." || gt == ".|." || gt == "." ||
                    gq == "." || gq == "" || gq + 0 < T[t] + 0) {
                    miss++
                    m[t, s] = 1
                } else {
                    m[t, s] = 0
                }
            }
            if (miss / ns <= GENO_PRE + 0) {
                kept[t]++
                for (s = 2; s <= NF; s++) if (m[t, s]) smiss[t, s]++
            }
        }
        # GQ histogram, split by hom-ref vs carrying an alt allele. The two
        # distributions are very different in a joint-genotyped WGS callset and
        # mixing them is what makes a single percentage misleading.
        for (s = 2; s <= NF; s++) {
            split($s, g, ":")
            if (g[2] == "." || g[2] == "") continue
            v = g[2] + 0; if (v > 99) v = 99
            if (g[1] == "0/0" || g[1] == "0|0") { href[int(v/10)]++; nhref++ }
            else if (g[1] ~ /[1-9]/)            { halt[int(v/10)]++; nhalt++ }
        }
        NS = NF - 1
    }
    END {
        printf "\n  GQ distribution by genotype class (sampled calls)\n"
        printf "  %-10s %12s %12s\n", "GQ bin", "hom-ref %", "alt-carrying %"
        for (b = 0; b <= 9; b++)
            printf "  %-10s %11.1f%% %11.1f%%\n", (b*10) "-" (b*10+9),
                   (nhref ? 100*href[b]/nhref : 0), (nhalt ? 100*halt[b]/nhalt : 0)
        printf "  n calls: hom-ref %d, alt-carrying %d\n", nhref, nhalt

        printf "\n  replay of the 01 filter chain (%d samples, %d variants sampled)\n", NS, nv
        printf "  %-4s %14s %10s %12s %12s %10s %12s\n",
               "T", "variants kept", "% kept", "mean s.miss", "max s.miss",
               "mind hits", "--mind needed"
        for (t = 1; t <= nt; t++) {
            tot = 0; mx = 0; hits = 0
            for (s = 2; s <= NS + 1; s++) {
                f = kept[t] ? smiss[t, s] / kept[t] : 0
                tot += f
                if (f > mx) mx = f
                if (f > MIND + 0) hits++
            }
            # 2x the worst observed sample, rounded up to two decimals. This is
            # a 1.6 Mb sample of easy euchromatin, so the genome-wide maximum
            # will be higher; a threshold sitting just above the observed max is
            # a coin flip, not a filter.
            need = int(mx * 200 + 0.999) / 100
            if (need < 0.02) need = 0.02
            printf "  %-4s %14d %9.1f%% %12.4f %12.4f %10d %12.2f\n",
                   T[t], kept[t], 100*kept[t]/nv, tot/NS, mx, hits, need
        }
        printf "\n  T           = --vcf-min-gq value (T=0 is no genotype filter)\n"
        printf "  variants    = surviving --geno %s after masking\n", GENO_PRE
        printf "  s.miss      = per-sample missingness among the survivors\n"
        printf "  mind hits   = samples --mind %s (MIND_KZ) would delete\n", MIND
        printf "  --mind needed = 2x the worst sample; set MIND_KZ at or above it\n"
    }' "${WORK_DIR}/gq_sample.txt" | tee "${OUT}/gq_profile.txt"

log ""
log "HOW TO READ IT. Compare each row against T=0."
log "  A few percent of variants lost and mind hits = 0: keep KZ_GQ_MIN there."
log "  Variant loss in the tens of percent, or any mind hits: that threshold"
log "  removes coverage, not errors. Drop to the largest T that stays clean."
log "  If even T=10 loses heavily, set KZ_GQ_MIN=0 and say so in Methods -"
log "  site-level PASS plus --geno/--mind is a defensible filter, and claiming"
log "  a genotype filter that silently deleted the panel is not."
log ""
log "Profile written to ${OUT}/gq_profile.txt (it belongs in the supplement)."
